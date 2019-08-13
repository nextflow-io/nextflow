package nextflow.script

import spock.lang.Specification

import nextflow.ast.NextflowDSL
import org.codehaus.groovy.control.CompilerConfiguration
import org.codehaus.groovy.control.customizers.ASTTransformationCustomizer
import test.MockSession
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ProcessDefTest extends Specification {

    def 'should define processes' () {

        given:
        def session = new MockSession()
        def binding = new ScriptBinding(session).setModule(true)
        def config = new CompilerConfiguration()
        config.setScriptBaseClass(BaseScript.class.name)
        config.addCompilationCustomizers( new ASTTransformationCustomizer(NextflowDSL))

        def SCRIPT = '''
            
            process foo {
              input: val data 
              output: val result
              exec:
                result = "$data mundo"
            }     
            
            process bar {
                input: val data 
                output: val result
                exec: 
                  result = data.toUpperCase()
            }   
            
        '''

        when:
        def script = (BaseScript)new GroovyShell(binding,config).parse(SCRIPT)

        then:
        true

    }

    def 'should clone a process with a new name - deprecated'() {

        given:
        def proc = new ProcessDef(Mock(BaseScript), 'foo', Mock(ProcessConfig), Mock(BodyDef))

        when:
        def copy = proc.cloneWithName('foo_alias')
        then:
        copy.getName() == 'foo_alias'

    }

    def 'should clone a process with a new name '() {

        given:
        def OWNER = Mock(BaseScript)
        def BODY = { -> null }
        def proc = new ProcessDef(OWNER, BODY, 'foo')

        when:
        def copy = proc.cloneWithName('foo_alias')
        then:
        copy.getName() == 'foo_alias'
        copy.getOwner() == OWNER
        copy.rawBody.class == BODY.class

    }
}
