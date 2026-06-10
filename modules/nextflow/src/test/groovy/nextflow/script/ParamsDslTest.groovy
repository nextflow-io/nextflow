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

import nextflow.Session
import nextflow.exception.AbortOperationException
import nextflow.exception.ScriptRuntimeException
import nextflow.file.FileHelper
import nextflow.script.types.Bag
import nextflow.script.types.Record
import spock.lang.Specification
import spock.lang.Unroll

import static test.ScriptHelper.*
/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class ParamsDslTest extends Specification {

    def 'should declare workflow params with CLI overrides'() {
        given:
        def inputFile = Files.createTempFile('test', '.csv')
        def cliParams = [input: inputFile.toString(), chunk_size: '3']
        def configParams = [outdir: 'results']

        when:
        def result = runScript(
            '''\
            params {
                input: Path
                chunk_size: Integer = 1
                save_intermeds: Boolean
            }

            workflow { params }
            ''',
            config: [params: configParams + cliParams],
            params: cliParams,
            configParams: configParams
        )
        then:
        result == [input: inputFile, chunk_size: 3, save_intermeds: false, outdir: 'results']

        cleanup:
        inputFile?.delete()
    }

    def 'should allow optional param'() {
        when:
        runScript(
            '''\
            params {
                input: Path?
            }

            workflow { params }
            '''
        )
        then:
        noExceptionThrown()
    }

    def 'should report error for missing required param'() {
        given:
        def cliParams = [:]
        def configParams = [outdir: 'results']

        when:
        runScript(
            '''\
            params {
                input: Path
                save_intermeds: Boolean
            }

            workflow { params }
            ''',
            params: cliParams,
            configParams: configParams
        )
        then:
        def e = thrown(ScriptRuntimeException)
        e.message == 'Parameter `input` is required but was not specified on the command line, params file, or config'
    }

    def 'should report error for invalid param'() {
        given:
        def cliParams = [inputs: './data']
        def configParams = [outdir: 'results']

        when:
        runScript(
            '''\
            params {
                input: Path
                save_intermeds: Boolean
            }

            workflow { params }
            ''',
            params: cliParams,
            configParams: configParams
        )
        then:
        def e = thrown(ScriptRuntimeException)
        e.message == 'Parameter `inputs` was specified on the command line or params file but is not declared in the script or config'
    }

    def 'should report error for invalid type'() {
        when:
        runScript(
            '''\
            params {
                save_intermeds: Boolean
            }

            workflow { params }
            ''',
            params: [save_intermeds: 42]
        )
        then:
        def e = thrown(ScriptRuntimeException)
        e.message == 'Parameter `save_intermeds` with type Boolean cannot be assigned to 42 [Integer]'
    }

    def 'should report error for missing input file'() {
        given:
        def cliParams = [input: 'input.csv']
        def configParams = [:]

        when:
        runScript(
            '''\
            params {
                input: Path
            }

            workflow { params }
            ''',
            params: cliParams,
            configParams: configParams
        )
        then:
        def e = thrown(AbortOperationException)
        e.message == "Input file 'input.csv' does not exist"
    }

    @Unroll
    def 'should validate float param with default value'() {
        when:
        runScript(
            """\
            params {
                factor: Float = ${DEF_VALUE}
            }

            workflow { params }
            """
        )
        then:
        noExceptionThrown()

        where:
        DEF_VALUE << [ '0.1f', '0.1d', '0.1g' ]
    }

    @Unroll
    def 'should validate integer param with default value'() {
        when:
        runScript(
            """\
            params {
                factor: Integer = ${DEF_VALUE}
            }

            workflow { params }
            """
        )
        then:
        noExceptionThrown()

        where:
        DEF_VALUE << [ '100i', '100l', '100g' ]
    }

    def 'should validate record collection param'() {
        given:
        def samples = [
            [id: 1, name: "sample1", value: 100],
            [id: 2, name: "sample2", value: 200],
            [id: 3, name: "sample3", value: 300]
        ]

        when:
        def result = runScript(
            '''\
            params {
                samples: List<Sample>
            }

            record Sample {
                id: Integer
                name: String
                value: Integer
            }

            workflow {
                params.samples
            }
            ''',
            params: [samples: samples]
        )
        then:
        result instanceof List
        result.size() == 3
        result[0] instanceof Record
        result[0].id == 1
        result[0].name == 'sample1'
        result[0].value == 100
        result[1].id == 2
        result[2].id == 3
    }

    def 'should report error for invalid record type'() {
        when:
        runScript(
            '''\
            params {
                items: List<Record>
            }

            workflow { params.items }
            ''',
            params: [items: [not: 'a list']]
        )
        then:
        def e = thrown(ScriptRuntimeException)
        e.message.contains('Parameter `items` with type List<Record> cannot be assigned to')
    }

    def 'should report error for missing record field'() {
        given:
        def samples = [
            [id: 1, name: "sample1"]
        ]

        when:
        runScript(
            '''\
            params {
                samples: List<Sample>
            }

            record Sample {
                id: Integer
                name: String
                value: Integer
            }

            workflow {
                params.samples
            }
            ''',
            params: [samples: samples]
        )

        then:
        def e = thrown(AbortOperationException)
        e.message == "Input record [id:1, name:sample1] is missing field 'value' required by record type 'Sample'"
    }

}
