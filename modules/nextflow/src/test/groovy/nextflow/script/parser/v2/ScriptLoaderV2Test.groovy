package nextflow.script.parser.v2

import java.nio.file.Files

import nextflow.Session
import nextflow.exception.ScriptCompilationException
import nextflow.script.BaseScript
import nextflow.script.ScriptMeta
import nextflow.script.WorkflowDef
import test.Dsl2Spec
import test.MockExecutorFactory
/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class ScriptLoaderV2Test extends Dsl2Spec {

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
        e.message == 'Script compilation failed'

        cleanup:
        file.parent.deleteDir()
    }

    def 'should register workflow definition' () {

        given:
        def session = new Session()
        def parser = new ScriptLoaderV2(session)

        def TEXT = '''

            workflow hello {
                take:
                foo ; bar

                main:
                def foobar = foo + bar

                emit:
                result = foobar
            }

            workflow {
                hello('foo', 'bar')
            }
            '''

        when:
        parser.parse(TEXT)
        parser.runScript()
        def meta = ScriptMeta.get(parser.getScript())

        then:
        meta.definitions.size() == 2
        meta.getWorkflow('hello').declaredInputs == ['foo', 'bar']
        meta.getWorkflow('hello').declaredOutputs == ['result']
    }

    def 'should register fully-qualified process names' () {

        given:
        def session = new Session()
        session.executorFactory = new MockExecutorFactory()
        def parser = new ScriptLoaderV2(session)

        def TEXT = '''

            process hello {
                'echo hello'
            }

            process bye {
                'echo bye'
            }

            workflow aloha {
                if( params.mode == 'hello' )
                    hello()
                if( params.mode == 'bye' )
                    bye()
            }

            workflow {
                aloha()
            }
            '''

        when:
        parser.parse(TEXT)
        parser.runScript()

        then:
        ScriptMeta.allProcessNames() == [ 'hello', 'bye', 'aloha:hello', 'aloha:bye' ] as Set
    }

    def 'should allow explicit `it` closure parameter' () {

        given:
        def session = new Session()
        def parser = new ScriptLoaderV2(session)

        def TEXT = '''
            workflow { 
                channel.of(1, 2, 3).view { it -> "${it}" }
            }
            '''

        when:
        parser.parse(TEXT)
        parser.runScript()

        then:
        noExceptionThrown()
    }

    def 'should eagerly evaluate GStrings' () {

        given:
        def session = new Session()
        def parser = new ScriptLoaderV2(session)

        def TEXT = '''
            workflow { 
                assert "${'hello'}" == 'hello'
                assert "${'hello'}" in ['hello']
            }
            '''

        when:
        parser.parse(TEXT)
        parser.runScript()

        then:
        noExceptionThrown()
    }

    def 'should lazily evaluate process inputs/outputs/directives' () {

        given:
        def session = new Session()
        session.executorFactory = new MockExecutorFactory()
        def parser = new ScriptLoaderV2(session)

        def TEXT = '''
            process HELLO {
                tag props.name

                input:
                val props

                output:
                val props.name

                script:
                "echo 'Hello ${props.name}'"
            }

            workflow {
                HELLO( [name: 'World'] )
            }
            '''

        when:
        parser.parse(TEXT)
        parser.runScript()

        then:
        noExceptionThrown()
    }

    def 'should strip unsupported type annotations' () {

        given:
        def session = new Session()
        def parser = new ScriptLoaderV2(session)

        def TEXT = '''
            // strip cast type
            ['1', '1.fastq', '2.fastq'] as Tuple<String,String,String>

            // strip type annotation in variable declaration
            def ch: Channel = channel.empty()
            '''

        when:
        parser.parse(TEXT)
        parser.runScript()

        then:
        noExceptionThrown()
    }

    def 'should register outputs and topic emissions' () {

        given:
        def session = new Session()
        def parser = new ScriptLoaderV2(session)

        def TEXT = '''

            nextflow.preview.types = true

            process hello {

                output:
                out1: Path = file("out1.txt")
                out2: Path = file("out2.txt")

                topic:
                file("version.txt") >> 'versions'

                script:
                """
                echo "version" > version.txt
                echo "result1" > out1.txt
                echo "result2" > out2.txt
                """
            }
            '''

        when:
        parser.setModule(true)
        parser.parse(TEXT)
        parser.runScript()
        def process = ScriptMeta.get(parser.getScript()).getProcess('hello')
        def outputs = process.getProcessConfig().getOutputs()

        then:
        outputs.getParams().size() == 2
        outputs.getTopics().size() == 1
        outputs.getFiles().size() == 3
    }

    def 'should support enums' () {

        given:
        def session = new Session()
        def parser = new ScriptLoaderV2(session)

        def TEXT = '''
            enum Day {
                MONDAY,
                TUESDAY,
                WEDNESDAY,
                THURSDAY,
                FRIDAY,
                SATURDAY,
                SUNDAY
            }

            workflow {
                Day.TUESDAY
            }
            '''

        when:
        parser.parse(TEXT)
        parser.runScript()

        then:
        parser.getResult().toString() == 'TUESDAY'
    }

}
