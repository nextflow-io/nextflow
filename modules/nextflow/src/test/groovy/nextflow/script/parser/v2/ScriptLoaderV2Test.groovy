package nextflow.script.parser.v2

import java.nio.file.Files

import nextflow.Session
import nextflow.exception.ScriptCompilationException
import nextflow.script.BaseScript
import nextflow.script.ScriptMeta
import nextflow.script.WorkflowDef
import spock.lang.Specification
/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class ScriptLoaderV2Test extends Specification {

    def 'should run a file script' () {

        given:
        def session = new Session()
        def parser = new ScriptLoaderV2(session)

        def file = Files.createTempDirectory('test').resolve('foo.nf')
        file.text = '''
        println "Hello world!"
        '''

        when:
        parser.parse(file)
        parser.runScript()
        then:
        parser.script instanceof BaseScript
        parser.binding.getScriptPath() == file
        parser.binding.getSession() == session
        and:
        def definitions = ScriptMeta.get(parser.script).getDefinitions()
        definitions.size() == 1
        definitions.first() instanceof WorkflowDef

        cleanup:
        file.parent.deleteDir()
    }

    def 'should run a text script' () {

        given:
        def session = new Session()
        def parser = new ScriptLoaderV2(session)

        def TEXT = '''
        println "Hello world!"
        '''

        when:
        parser.runScript(TEXT)
        then:
        parser.script instanceof BaseScript
        parser.binding.getScriptPath() == null
        parser.binding.getSession() == session
        and:
        def definitions = ScriptMeta.get(parser.script).getDefinitions()
        definitions.size() == 1
        definitions.first() instanceof WorkflowDef
    }

    def 'should report compilation errors' () {
        given:
        def session = new Session()
        def parser = new ScriptLoaderV2(session)

        def file = Files.createTempDirectory('test').resolve('foo.nf')
        file.text = '''
        def foo)( { )
        '''

        when:
        parser.parse(file)
        parser.runScript()
        then:
        def e = thrown(ScriptCompilationException)
        e.message.contains("Unexpected input: '{'")
        e.message.contains('foo.nf')

        cleanup:
        file.parent.deleteDir()
    }
    
}
