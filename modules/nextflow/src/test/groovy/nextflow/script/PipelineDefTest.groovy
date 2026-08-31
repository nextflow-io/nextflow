/*
 * Copyright 2013-2026, Seqera Labs
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package nextflow.script

import java.nio.file.Files
import java.nio.file.Path

import nextflow.exception.ScriptCompilationException
import nextflow.exception.ScriptRuntimeException
import spock.lang.Timeout
import test.Dsl2Spec

import static test.ScriptHelper.*

/**
 * Tests for {@link PipelineDef} -- a pipeline included as a named workflow.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Timeout(30)
class PipelineDefTest extends Dsl2Spec {

    private Path folder

    def setup() {
        folder = Files.createTempDirectory('test')
    }

    def cleanup() {
        folder?.deleteDir()
    }

    private static final String GREET_PIPELINE = """
            nextflow.enable.types = true

            params {
                names: Channel<String>
                greeting: String = 'Hello'
                subdir: String = 'greetings'
            }

            workflow {
                main:
                messages = params.names.map { name -> "\${params.greeting}, \${name}!" }

                publish:
                messages = messages
            }

            output {
                messages: Channel<String> { path params.subdir }
            }
            """

    private Path pipeline(String text) {
        folder.resolve('greet.nf').text = """
            nextflow.enable.types = true

            params {
                names: Channel<String>
                greeting: String = 'Hello'
            }

            workflow {
                main:
                messages = params.names.map { name -> "\${params.greeting}, \${name}!" }

                publish:
                messages = messages
            }

            output {
                messages: Channel<String> {}
            }
            """

        folder.resolve('main.nf').text = """
            nextflow.enable.types = true

            include {
                params as GreetParams ;
                workflow as GREET ;
                output as GreetOutput
            } from './greet.nf'

            ${text}
            """

        return folder.resolve('main.nf')
    }

    def 'should call an included pipeline like a named workflow' () {
        given:
        def script = pipeline('''
            workflow {
                main:
                greet = GREET( names: channel.of('World', 'Nextflow') )
                greet.messages
            }
            ''')

        when:
        def result = runScript(script)
        then:
        result.val == 'Hello, World!'
        result.val == 'Hello, Nextflow!'
    }

    def 'should override a param default' () {
        given:
        def script = pipeline('''
            workflow {
                main:
                greet = GREET( names: channel.of('World'), greeting: 'Hola' )
                greet.messages
            }
            ''')

        when:
        def result = runScript(script)
        then:
        result.val == 'Hola, World!'
    }

    def 'should accept the params as a record' () {
        given:
        def script = pipeline('''
            workflow {
                main:
                opts = record( greeting: 'Ciao' )
                greet = GREET( opts + record( names: channel.of('World') ) )
                greet.messages
            }
            ''')

        when:
        def result = runScript(script)
        then:
        result.val == 'Ciao, World!'
    }

    def 'should fail when a required param is not provided' () {
        given:
        def script = pipeline('''
            workflow {
                main:
                GREET( greeting: 'Hola' )
            }
            ''')

        when:
        runScript(script)
        then:
        def e = thrown(ScriptRuntimeException)
        e.message == 'Parameter `names` of pipeline `GREET` is required but no value was provided'
    }

    def 'should fail when an undeclared param is provided' () {
        given:
        def script = pipeline('''
            workflow {
                main:
                GREET( names: channel.of('World'), foo: 'bar' )
            }
            ''')

        when:
        runScript(script)
        then:
        def e = thrown(ScriptRuntimeException)
        e.message == 'Pipeline `GREET` does not declare a parameter named `foo`'
    }

    def 'should call an included pipeline with its params record' () {
        given:
        def script = pipeline("""
            params {
                greet: GreetParams
            }

            workflow {
                main:
                greet = GREET( params.greet + record(names: channel.of('World')) )
                greet.messages
            }
            """)

        when:
        def result = runScript([params: [greet: [greeting: 'Hola']]], script)
        then:
        result.val == 'Hola, World!'
    }

    def 'should apply the pipeline default for a param omitted from its params record' () {
        given:
        def script = pipeline("""
            params {
                greet: GreetParams
            }

            workflow {
                main:
                greet = GREET( params.greet + record(names: channel.of('World')) )
                greet.messages
            }
            """)

        when:
        def result = runScript([params: [greet: [:]]], script)
        then:
        result.val == 'Hello, World!'
    }

    def 'should publish the outputs of an included pipeline from its output record' () {
        given:
        def script = pipeline("""
            workflow {
                main:
                greet = GREET( names: channel.of('World') )

                publish:
                greet = greet
            }

            output {
                greet: GreetOutput {}
            }
            """)

        when:
        // the output block declares one output for each output of the included
        // pipeline, and the record of channels is published field by field --
        // a mismatch on either side is reported when the output block is applied
        runScript([config: [outputDir: folder.resolve('results').toString()]], script)
        then:
        noExceptionThrown()
    }

    def 'should scope the processes of an included pipeline by its name' () {
        given:
        folder.resolve('greet.nf').text = '''
            nextflow.enable.types = true

            params {
                names: Channel<String>
            }

            workflow {
                main:
                messages = FOO(params.names)

                publish:
                messages = messages
            }

            output {
                messages: Channel<String> {}
            }

            process FOO {
                input:
                name: String

                output:
                message: String

                exec:
                message = "${task.ext.greeting}, ${name}!"
            }
            '''

        folder.resolve('main.nf').text = '''
            nextflow.enable.types = true

            include { workflow as GREET } from './greet.nf'

            workflow {
                main:
                greet = GREET( names: channel.of('World') )
                greet.messages
            }
            '''

        when:
        def result = runScript([config: [process: ['withName:GREET:FOO': [ext: [greeting: 'Bonjour']]]]], folder.resolve('main.nf'))
        then:
        result.val == 'Bonjour, World!'
    }

    def 'should publish the outputs of the calling pipeline only' () {
        given:
        def script = pipeline('''
            workflow {
                main:
                greet = GREET( names: channel.of('World') )

                publish:
                out = greet.messages
            }

            output {
                out: Channel<String> {}
            }
            ''')

        when:
        runScript(script)
        then:
        noExceptionThrown()
    }

    def 'should fail when an undeclared param is provided as null' () {
        given:
        def script = pipeline('''
            workflow {
                main:
                GREET( names: channel.of('World'), foo: null )
            }
            ''')

        when:
        runScript(script)
        then:
        def e = thrown(ScriptRuntimeException)
        e.message == 'Pipeline `GREET` does not declare a parameter named `foo`'
    }

    def 'should not expose the params of the calling pipeline to an included pipeline' () {
        given:
        folder.resolve('greet.nf').text = '''
            nextflow.enable.types = true

            params {
                names: Channel<String>
            }

            workflow {
                main:
                messages = params.names.map { name -> "${params.secret}, ${name}!" }

                publish:
                messages = messages
            }

            output {
                messages: Channel<String> {}
            }
            '''

        folder.resolve('main.nf').text = '''
            nextflow.enable.types = true

            include { workflow as GREET } from './greet.nf'

            params {
                secret: String = 'LEAKED'
            }

            workflow {
                main:
                greet = GREET( names: channel.of('World') )
                greet.messages
            }
            '''

        when:
        def result = runScript(folder.resolve('main.nf'))
        then:
        // a param that the included pipeline does not declare is not visible
        // to it, even though the calling pipeline declares it
        result.val == 'null, World!'
    }

    def 'should not require a params record to be provided at launch' () {
        given:
        def script = pipeline('''
            params {
                greet: GreetParams
            }

            workflow {
                main:
                greet = GREET( params.greet + record(names: channel.of('World')) )
                greet.messages
            }
            ''')

        when:
        // a partial record type defaults to an empty record, so a pipeline
        // whose params are supplied by the calling pipeline can be launched
        // without providing any of them
        def result = runScript(script)
        then:
        result.val == 'Hello, World!'
    }

    def 'should merge a nested command-line param with the config value' () {
        given:
        def script = pipeline('''
            params {
                greet: GreetParams
            }

            workflow {
                main:
                greet = GREET( params.greet + record(names: channel.of('World')) )
                greet.messages
            }
            ''')

        when:
        // the command line overrides only the field that it names
        def result = runScript([params: [greet: [greeting: 'Hola']], configParams: [greet: [greeting: 'Ciao']]], script)
        then:
        result.val == 'Hola, World!'
    }

    def 'should support multiple aliases of the same pipeline' () {
        given:
        def script = pipeline('''
            include { workflow as GREET_AGAIN } from './greet.nf'

            workflow {
                main:
                a = GREET( names: channel.of('World'), greeting: 'Hola' )
                b = GREET_AGAIN( names: channel.of('Nextflow'), greeting: 'Ciao' )
                a.messages.mix(b.messages)
            }
            ''')

        when:
        def result = runScript(script)
        then:
        // each alias resolves its own params, because the entry workflow
        // receives them as an input instead of reading global state
        [result.val, result.val].sort() == ['Ciao, Nextflow!', 'Hola, World!']
    }

    def 'should support calling the same alias more than once when the pipeline has no process' () {
        given:
        def script = pipeline('''
            workflow {
                main:
                a = GREET( names: channel.of('World'), greeting: 'Hola' )
                b = GREET( names: channel.of('Nextflow'), greeting: 'Ciao' )
                a.messages.mix(b.messages)
            }
            ''')

        when:
        def result = runScript(script)
        then:
        // a pipeline that contains a process can be called only once per
        // alias, like a named workflow -- see the spec below
        [result.val, result.val].sort() == ['Ciao, Nextflow!', 'Hola, World!']
    }

    def 'should reject calling the same alias more than once when the pipeline has a process' () {
        given:
        folder.resolve('greet.nf').text = '''
            nextflow.enable.types = true

            params {
                greeting: String = 'Hello'
            }

            process SAY {
                input:
                message: String
                output:
                stdout()
                script:
                "echo '${message}'"
            }

            workflow {
                main:
                said = SAY( params.greeting )

                publish:
                said = said
            }

            output {
                said: Channel<String> {}
            }
            '''
        folder.resolve('main.nf').text = '''
            nextflow.enable.types = true

            include { workflow as GREET } from './greet.nf'

            workflow {
                main:
                GREET( greeting: 'Hola' )
                GREET( greeting: 'Ciao' )
            }
            '''

        when:
        runScript(folder.resolve('main.nf'))
        then:
        def e = thrown(Exception)
        e.message.contains("has been already used")
    }

    def 'should reject a pipeline that refers to params outside its entry workflow' () {
        given:
        folder.resolve('greet.nf').text = '''
            nextflow.enable.types = true

            params {
                names: Channel<String>
                greeting: String = 'Hello'
            }

            workflow {
                main:
                messages = FOO(params.names)

                publish:
                messages = messages
            }

            output {
                messages: Channel<String> {}
            }

            process FOO {
                input:
                name: String

                output:
                message: String

                exec:
                message = "${params.greeting}, ${name}!"
            }
            '''

        folder.resolve('main.nf').text = '''
            nextflow.enable.types = true

            include { workflow as GREET } from './greet.nf'

            workflow {
                main:
                GREET( names: channel.of('World') )
            }
            '''

        when:
        runScript(folder.resolve('main.nf'))
        then:
        def e = thrown(ScriptCompilationException)
        e.cause.message.contains('An included pipeline cannot refer to `params` outside of its entry workflow and output block')
    }

    def 'should resolve the output directives of a pipeline against its params' () {
        when:
        def script = loadScript(module: true, GREET_PIPELINE)
        then:
        // the output block is evaluated per pipeline execution, so its
        // directives can refer to the params of that execution
        script.getOutputDeclarations(params(subdir: 'hola')).messages.path == 'hola'
        script.getOutputDeclarations(params(subdir: 'ciao')).messages.path == 'ciao'
    }

    private static ScriptBinding.ParamsMap params(Map values) {
        return new ScriptBinding.ParamsMap(values)
    }

    def 'should fail when an included pipeline is not aliased' () {
        given:
        def script = pipeline('workflow { }')
        folder.resolve('main.nf').text = folder.resolve('main.nf').text
            .replace('workflow as GREET', 'workflow')

        when:
        runScript(script)
        then:
        def e = thrown(ScriptCompilationException)
        e.cause.message.contains('An included pipeline must be aliased')
    }


    def 'should support including the same params block from two scripts' () {
        given:
        def script = pipeline('''
            include { params as MidParams ; workflow as MID } from './mid.nf'

            workflow {
                main:
                MID( names: channel.of('World') ).messages
            }
            ''')
        folder.resolve('mid.nf').text = '''
            nextflow.enable.types = true

            include { params as GreetParams ; workflow as GREET } from './greet.nf'

            params {
                names: Channel<String>
            }

            workflow {
                main:
                messages = GREET( names: params.names ).messages

                publish:
                messages = messages
            }

            output {
                messages: Channel<String> {}
            }
            '''

        when:
        // the record type of an included params block is qualified by the
        // including script, so two scripts can include the same one
        def result = runScript(script)
        then:
        result.val == 'Hello, World!'
    }

    def 'should execute a pipeline module only once when it is included by two scripts' () {
        given:
        def script = pipeline('''
            include { workflow as MID } from './mid.nf'

            workflow {
                main:
                a = GREET( names: channel.of('World'), greeting: 'Hola' )
                b = MID( names: channel.of('Nextflow') )
                a.messages.mix(b.messages)
            }
            ''')
        folder.resolve('mid.nf').text = '''
            nextflow.enable.types = true

            include { workflow as GREET_INNER } from './greet.nf'

            params {
                names: Channel<String>
            }

            workflow {
                main:
                messages = GREET_INNER( names: params.names, greeting: 'Ciao' ).messages

                publish:
                messages = messages
            }

            output {
                messages: Channel<String> {}
            }
            '''

        when:
        def result = runScript(script)
        then:
        [result.val, result.val].sort() == ['Ciao, Nextflow!', 'Hola, World!']
    }

    def 'should reject a pipeline whose module refers to params' () {
        given:
        folder.resolve('greet.nf').text = '''
            nextflow.enable.types = true

            include { SUB } from './sub.nf'

            params {
                names: Channel<String>
            }

            workflow {
                main:
                messages = SUB( params.names )

                publish:
                messages = messages
            }

            output {
                messages: Channel<String> {}
            }
            '''
        folder.resolve('sub.nf').text = '''
            nextflow.enable.types = true

            workflow SUB {
                take:
                names: Channel<String>

                main:
                messages = names.map { name -> "${params.greeting}, ${name}!" }

                emit:
                messages
            }
            '''
        folder.resolve('main.nf').text = '''
            nextflow.enable.types = true

            include { workflow as GREET } from './greet.nf'

            workflow {
                main:
                GREET( names: channel.of('World') )
            }
            '''

        when:
        // a `params` reference in a module of the pipeline would resolve
        // against the calling pipeline, so it is rejected as well
        runScript(folder.resolve('main.nf'))
        then:
        def e = thrown(ScriptCompilationException)
        e.cause.message.contains('An included pipeline cannot refer to `params` outside of its entry workflow and output block')
    }

    def 'should allow a pipeline that binds a local variable named params' () {
        given:
        folder.resolve('greet.nf').text = '''
            nextflow.enable.types = true

            params {
                names: Channel<String>
            }

            def render(params: Map) -> String {
                return params.toString()
            }

            workflow {
                main:
                messages = params.names.map { name -> render([name: name]) }

                publish:
                messages = messages
            }

            output {
                messages: Channel<String> {}
            }
            '''
        folder.resolve('main.nf').text = '''
            nextflow.enable.types = true

            include { workflow as GREET } from './greet.nf'

            workflow {
                main:
                GREET( names: channel.of('World') ).messages
            }
            '''

        when:
        // `params` is shadowed by the function input, so it is not a
        // reference to the params of the pipeline
        def result = runScript(folder.resolve('main.nf'))
        then:
        result.val == '[name:World]'
    }

    def 'should include a definition named output' () {
        given:
        folder.resolve('lib.nf').text = '''
            nextflow.enable.types = true

            def output(value: String) -> String {
                return value.toUpperCase()
            }
            '''
        folder.resolve('main.nf').text = '''
            nextflow.enable.types = true

            include { output } from './lib.nf'

            workflow {
                main:
                channel.of('World').map { name -> output(name) }
            }
            '''

        when:
        // a definition that happens to be named `output` is not the output
        // block of a pipeline, so it needs no alias
        def result = runScript(folder.resolve('main.nf'))
        then:
        result.val == 'WORLD'
    }

}
