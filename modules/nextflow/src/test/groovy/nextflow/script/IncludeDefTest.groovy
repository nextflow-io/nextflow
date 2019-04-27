package nextflow.script

import spock.lang.Specification
import spock.lang.Unroll

import java.nio.file.NoSuchFileException
import java.nio.file.Path

import nextflow.ast.NextflowDSL
import org.codehaus.groovy.control.CompilerConfiguration
import org.codehaus.groovy.control.customizers.ASTTransformationCustomizer
import test.TestHelper
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class IncludeDefTest extends Specification {

    @Unroll
    def 'should resolve module path #MODULE' () {

        given:
        def script = '/some/path/main.nf' as Path
        def include = Spy(IncludeDef)
        include.getOwnerPath() >> script

        expect:
        include.resolveModulePath('/abs/foo.nf') == '/abs/foo.nf' as Path
        include.resolveModulePath('module.nf') == '/some/path/module.nf' as Path
        include.resolveModulePath('foo/bar.nf') == '/some/path/foo/bar.nf' as Path

        when:
        include.resolveModulePath('http://foo.com/bar')
        then:
        thrown(IllegalArgumentException)

    }

    def 'should resolve real module path' () {

        given:
        def folder = TestHelper.createInMemTempDir()
        def script = folder.resolve('main.nf'); script.text = 'echo ciao'
        def module = folder.resolve('mod-x.nf'); module.text = 'blah blah'

        def include = Spy(IncludeDef)
        include.getOwnerPath() >> script
        
        when:
        def result = include.realModulePath( 'mod-x.nf')
        then:
        result == module

        when:
        result = include.realModulePath('mod-x')
        then:
        result == module

        when:
        include.realModulePath('xyz')
        then:
        thrown(NoSuchFileException)

    }

    // ==== DSL tests === 

    static class TestInclude extends IncludeDef {
        boolean loadInvoked

        TestInclude(IncludeDef include) {
            this.path = include.path
            this.name = include.name
            this.alias = include.alias
        }

        @Override
        void load() {
            loadInvoked = true
        }
    }

    static abstract class TestScript extends BaseScript {

        List<TestInclude> includes = []

        @Override
        protected IncludeDef include(IncludeDef include) {
            def inc = new TestInclude(include)
            includes << inc
            return inc
        }

        @Override
        Object run() {
            return runScript()
        }
    }

    def 'should add includes' () {
        given:
        def config = new CompilerConfiguration()
        config.setScriptBaseClass(TestScript.class.name)
        config.addCompilationCustomizers( new ASTTransformationCustomizer(NextflowDSL))

        when:
        def script = (TestScript)new GroovyShell(config).parse(INCLUDE)
        script.run()
        then:
        script.includes[0].path == PATH
        script.includes[0].name == NAME
        script.includes[0].alias == ALIAS
        script.includes[0].params == PARAMS
        script.includes[0].loadInvoked


        where:
        INCLUDE                                         | PATH          | NAME      | ALIAS     | PARAMS
        "include 'some/path'"                           | 'some/path'   | null      | null      | [:]
        "include ALPHA from 'modules/path'"             | 'modules/path'| 'ALPHA'   | null      | [:]
        "include ALPHA as BRAVO from 'modules/x'"       | 'modules/x'   | 'ALPHA'   | 'BRAVO'   | [:]
        "include 'modules/1' params(a:1, b:2)"          | 'modules/1'   | null      | null      | [a:1, b:2]
        "include DELTA from 'abc' params(x:1)"          | 'abc'         | 'DELTA'   | null      | [x:1]
        "include GAMMA as FOO from 'm1' params(p:2)"    | 'm1'          | 'GAMMA'   | 'FOO'     | [p:2]

    }


}
