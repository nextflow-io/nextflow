package nextflow.script

import spock.lang.Specification

import java.nio.file.Paths

import nextflow.Session
import test.TestHelper
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ScriptParserTest extends Specification {

    def 'should run a file script' () {

        given:
        def session = new Session()
        def parser = new ScriptParser(session)
        def binding = new ScriptBinding(params:[foo:'Hello'])

        def file = TestHelper.createInMemTempFile('foo.nf')
        file.text = '''
        bar = "$params.foo world!"
        '''

        when:
        parser.setBinding(binding)
        parser.runScript(file)
        then:
        parser.script instanceof BaseScript
        parser.result == 'Hello world!'
        parser.result == binding.getVariable('bar')
        parser.binding.scriptPath == file
        parser.binding.module == false
        parser.binding.session == session
        !session.binding.hasVariable('bar')
    }

    def 'should run a text script' () {

        given:
        def session = new Session()
        def parser = new ScriptParser(session)
        def binding = new ScriptBinding(params:[foo:'Hello'])

        def TEXT = '''
        bar = "$params.foo world!"
        '''

        when:
        parser.setBinding(binding)
        parser.runScript(TEXT)
        then:
        parser.script instanceof BaseScript
        parser.result == 'Hello world!'
        parser.result == binding.getVariable('bar')
        parser.binding.scriptPath == null
        parser.binding.module == false
        parser.binding.session == session
        !session.binding.hasVariable('bar')
    }

    def 'should run a script with session binding' () {

        given:
        def session = new Session()
        def parser = new ScriptParser(session)
        session.binding.setVariable('foo', 'Hello')

        def TEXT = '''
        bar = "$foo world!"
        '''

        when:
        parser.runScript(TEXT)
        then:
        parser.script instanceof BaseScript
        parser.result == 'Hello world!'
        parser.binding.scriptPath == null
        parser.binding.module == false
        parser.binding.session == session
        session.binding.getVariable('foo') == 'Hello'
        session.binding.getVariable('bar') == 'Hello world!'
    }

    def 'should normalise script name'() {

        given:
        def parser = new ScriptParser(Mock(Session))

        expect:
        parser.computeClassName(Paths.get(SCRIPT)) == EXPECTED

        where:
        SCRIPT              | EXPECTED
        'foo.nf'            | 'Script_foo'
        'foo-bar-baz.nf'    | 'Script_foo_bar_baz'
        '123-fo0'           | 'Script_23_fo0'
        '--a  b  c'         | 'Script_a_b_c'
    }

    def 'should normalise script text' () {

        given:
        def parser = new ScriptParser(Mock(Session))

        when:
        def result = parser.computeClassName('process foo { etc } ')
        then:
        result == 'Script_01af1441'
    }


}
