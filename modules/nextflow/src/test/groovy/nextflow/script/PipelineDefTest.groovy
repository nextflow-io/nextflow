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

import nextflow.exception.ScriptRuntimeException
import spock.lang.Timeout
import spock.lang.Unroll
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

    private static final String GREET = '''
        params {
            names: Channel<String>
            greeting: String = 'Hello'
        }

        workflow {
            main:
            messages = params.names.map { name -> "${params.greeting}, ${name}!" }

            publish:
            messages = messages
        }

        output {
            messages: Channel<String> {}
        }
        '''

    private static final String COUNT = '''
        params {
            samples: Channel<Sample>
            factor: Value<Integer> = 1
        }

        record Sample {
            id: String
            count: Integer
        }

        workflow {
            main:
            totals = params.samples
                .map { s -> s.count }
                .collect()
                .combine(params.factor)
                .map { counts, factor -> counts.sum() * factor }

            publish:
            total = totals
        }

        output {
            total: Channel<Integer> {}
        }
        '''

    /**
     * Write typed scripts to the test folder and return the path of `main.nf`.
     */
    private Path write(Map<String,String> files) {
        files.each { name, text ->
            folder.resolve(name).text = 'nextflow.enable.types = true\n' + text.stripIndent()
        }
        return folder.resolve('main.nf')
    }

    /**
     * Write the greet pipeline and a main script that includes it.
     */
    private Path pipeline(String text) {
        return write([
            'greet.nf': GREET,
            'main.nf': '''
                include {
                    params as GreetParams ;
                    workflow as GREET ;
                    output as GreetOutput
                } from './greet.nf'
                ''' + text
        ])
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

    @Unroll
    def 'should fail on an invalid pipeline call: #CALL' () {
        given:
        def script = pipeline("""
            workflow {
                main:
                ${CALL}
            }
            """)

        when:
        runScript(script)
        then:
        def e = thrown(ScriptRuntimeException)
        e.message == ERROR

        where:
        CALL                                                                    | ERROR
        "GREET( greeting: 'Hola' )"                                             | 'Parameter `names` of pipeline `GREET` is required but no value was provided'
        "GREET( names: channel.of('World'), greeting: null )"                   | 'Parameter `greeting` of pipeline `GREET` is required but no value was provided'
        "GREET( names: channel.of('World'), foo: 'bar' )"                       | 'Pipeline `GREET` does not declare a parameter named `foo`'
        "GREET( names: channel.of('World'), foo: null )"                        | 'Pipeline `GREET` does not declare a parameter named `foo`'
        "GREET( names: channel.of('World'), greeting: channel.value('Hola') )"  | 'Parameter `greeting` of pipeline `GREET` with type String cannot be assigned to a dataflow value -- declare the param as a Channel or Value to accept it'
        "GREET( names: channel.value('World') )"                                | 'Parameter `names` of pipeline `GREET` with type Channel<String> cannot be assigned to a Value'
        "GREET( names: 'World' )"                                               | 'Parameter `names` of pipeline `GREET` with type Channel<String> cannot be assigned to World [String]'
    }

    @Unroll
    def 'should resolve an included params record from the command line and config' () {
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
        def result = runScript([params: CLI, configParams: CONFIG], script)
        then:
        result.val == RESULT

        where:
        CLI                          | CONFIG                        | RESULT
        [greet: [greeting: 'Hola']]  | [:]                           | 'Hola, World!'
        [:]                          | [greet: [greeting: 'Ciao']]   | 'Ciao, World!'
        // the command line overrides only the field that it names
        [greet: [greeting: 'Hola']]  | [greet: [greeting: 'Ciao']]   | 'Hola, World!'
        // the pipeline applies its own default for an omitted param
        [greet: [:]]                 | [:]                           | 'Hello, World!'
        // a params record with no required fields is not required at launch
        [:]                          | [:]                           | 'Hello, World!'
    }

    def 'should include the output block of a pipeline as a record type' () {
        given:
        def script = pipeline('''
            workflow SHOUT {
                take:
                greet: GreetOutput

                main:
                messages = greet.messages.map { message -> message.toUpperCase() }

                emit:
                messages: Channel<String> = messages
            }

            workflow {
                main:
                SHOUT(GREET( names: channel.of('World') ))
            }
            ''')

        when:
        // the output block is imported as a record type of the pipeline
        // outputs -- it declares no outputs of its own
        def result = runScript(script)
        then:
        result.val == 'HELLO, WORLD!'
    }

    @Unroll
    def 'should accept dataflow values for Channel and Value params: #CALL' () {
        given:
        def script = write([
            'count.nf': COUNT,
            'main.nf': """
                include { workflow as COUNT } from './count.nf'

                workflow {
                    main:
                    samples = channel.of( record(id: 'a', count: 1), record(id: 'b', count: 2) )
                    ${CALL}.total
                }
                """
        ])

        when:
        def result = runScript(script)
        then:
        result.val == RESULT

        where:
        CALL                                                    | RESULT
        "COUNT( samples: samples, factor: channel.value(2) )"   | 6
        // the default of a Value param is wrapped in a dataflow value
        "COUNT( samples: samples )"                             | 3
    }

    def 'should fail when a channel is provided for a Value param' () {
        given:
        def script = write([
            'count.nf': COUNT,
            'main.nf': '''
                include { workflow as COUNT } from './count.nf'

                workflow {
                    main:
                    samples = channel.of( record(id: 'a', count: 1) )
                    COUNT( samples: samples, factor: channel.of(1, 2) )
                }
                '''
        ])

        when:
        runScript(script)
        then:
        def e = thrown(ScriptRuntimeException)
        e.message == 'Parameter `factor` of pipeline `COUNT` with type Value<Integer> cannot be assigned to a Channel'
    }

    def 'should load the dataflow params of an included pipeline from the command line' () {
        given:
        folder.resolve('samples.csv').text = 'id,count\na,1\nb,2\n'
        def script = write([
            'count.nf': COUNT,
            'main.nf': '''
                include { params as CountParams ; workflow as COUNT } from './count.nf'

                params {
                    count: CountParams
                }

                workflow {
                    main:
                    COUNT( params.count ).total
                }
                '''
        ])
        def cliParams = [count: [samples: folder.resolve('samples.csv').toString(), factor: '5']]

        when:
        def result = runScript([params: cliParams], script)
        then:
        result.val == 15
    }

    def 'should call a typed pipeline from an untyped script' () {
        given:
        write(['greet.nf': GREET])
        def script = folder.resolve('main.nf')
        script.text = '''
            include { workflow as GREET } from './greet.nf'

            workflow {
                GREET( names: channel.of('World') ).messages
            }
            '''

        when:
        def result = runScript(script)
        then:
        result.val == 'Hello, World!'
    }

    def 'should call an untyped pipeline from a typed script' () {
        given:
        folder.resolve('greet.nf').text = '''
            params {
                name: String
                greeting: String = 'Hello'
            }

            workflow {
                main:
                messages = channel.of("${params.greeting}, ${params.name}!")

                publish:
                messages = messages
            }

            output {
                messages {}
            }
            '''
        def script = write([
            'main.nf': '''
                include { workflow as GREET } from './greet.nf'

                workflow {
                    main:
                    GREET( name: 'World' ).messages
                }
                '''
        ])

        when:
        def result = runScript(script)
        then:
        result.val == 'Hello, World!'
    }

    def 'should scope the processes of an included pipeline by its name' () {
        given:
        def script = write([
            'greet.nf': '''
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
                ''',
            'main.nf': '''
                include { workflow as GREET } from './greet.nf'

                workflow {
                    main:
                    greet = GREET( names: channel.of('World') )
                    greet.messages
                }
                '''
        ])

        when:
        def result = runScript([config: [process: ['withName:GREET:FOO': [ext: [greeting: 'Bonjour']]]]], script)
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

    def 'should not expose the params of the calling pipeline to an included pipeline' () {
        given:
        def script = write([
            'greet.nf': '''
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
                ''',
            'main.nf': '''
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
        ])

        when:
        def result = runScript(script)
        then:
        // a param that the included pipeline does not declare is not visible
        // to its entry workflow, even though the calling pipeline declares it
        result.val == 'null, World!'
    }

    def 'should call an included pipeline without a params block' () {
        given:
        def script = write([
            'hello.nf': '''
                workflow {
                    main:
                    messages = channel.of('Hello')

                    publish:
                    messages = messages
                }

                output {
                    messages: Channel<String> {}
                }
                ''',
            'main.nf': '''
                include { workflow as HELLO } from './hello.nf'

                params {
                    greeting: String = 'Hola'
                }

                workflow {
                    main:
                    HELLO().messages
                }
                '''
        ])

        when:
        def result = runScript(script)
        then:
        result.val == 'Hello'
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
        def script = write([
            'greet.nf': '''
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
                ''',
            'main.nf': '''
                include { workflow as GREET } from './greet.nf'

                workflow {
                    main:
                    GREET( greeting: 'Hola' )
                    GREET( greeting: 'Ciao' )
                }
                '''
        ])

        when:
        runScript(script)
        then:
        def e = thrown(Exception)
        e.message.contains("Process 'GREET:SAY' was called twice")
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
        write([
            'mid.nf': '''
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
        ])

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
        write([
            'mid.nf': '''
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
        ])

        when:
        def result = runScript(script)
        then:
        [result.val, result.val].sort() == ['Ciao, Nextflow!', 'Hola, World!']
    }

    def 'should include a definition named output' () {
        given:
        def script = write([
            'lib.nf': '''
                def output(value: String) -> String {
                    return value.toUpperCase()
                }
                ''',
            'main.nf': '''
                include { output } from './lib.nf'

                workflow {
                    main:
                    channel.of('World').map { name -> output(name) }
                }
                '''
        ])

        when:
        // a definition that happens to be named `output` is not the output
        // block of a pipeline, so it needs no alias
        def result = runScript(script)
        then:
        result.val == 'WORLD'
    }

}
