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
        given:
        def cliParams = [save_intermeds: 42]
        def configParams = [:]

        when:
        runScript(
            '''\
            params {
                save_intermeds: Boolean
            }

            workflow { params }
            ''',
            params: cliParams,
            configParams: configParams
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

    def 'should load collection param from CSV file'() {
        given:
        def csvFile = Files.createTempFile('test', '.csv')
        csvFile.text = '''\
            id,name,value
            1,sample1,100
            2,sample2,200
            3,sample3,300
            '''.stripIndent()
        def cliParams = [samples: csvFile.toString()]

        when:
        def samples = runScript(
            '''\
            params {
                samples: List<Record>
            }

            workflow { params.samples }
            ''',
            params: cliParams,
            configParams: [:]
        )

        then:
        samples instanceof List
        samples.size() == 3
        samples[0].id == '1'
        samples[0].name == 'sample1'
        samples[0].value == '100'
        samples[1].id == '2'
        samples[2].id == '3'

        cleanup:
        csvFile?.delete()
    }

    def 'should load collection param from JSON file'() {
        given:
        def jsonFile = Files.createTempFile('test', '.json')
        jsonFile.text = '''\
            [
              {"id": 1, "name": "sample1", "value": 100},
              {"id": 2, "name": "sample2", "value": 200},
              {"id": 3, "name": "sample3", "value": 300}
            ]
            '''.stripIndent()
        def cliParams = [
            samplesList: jsonFile.toString(),
            samplesBag: jsonFile.toString(),
            samplesSet: jsonFile.toString()
        ]

        when:
        def params = runScript(
            '''\
            params {
                samplesList: List<Record>
                samplesBag: Bag<Record>
                samplesSet: Set<Record>
            }

            workflow { params }
            ''',
            params: cliParams,
            configParams: [:]
        )

        then:
        def samplesList = params.samplesList
        samplesList instanceof List
        samplesList.size() == 3
        samplesList[0].id == 1
        samplesList[0].name == 'sample1'
        samplesList[0].value == 100
        samplesList[1].id == 2
        samplesList[2].id == 3

        def samplesBag = params.samplesBag
        samplesBag instanceof Bag
        samplesBag.size() == 3

        def samplesSet = params.samplesSet
        samplesSet instanceof Set
        samplesSet.size() == 3

        cleanup:
        jsonFile?.delete()
    }

    def 'should load collection param from YAML file'() {
        given:
        def yamlFile = Files.createTempFile('test', '.yml')
        yamlFile.text = '''\
            - id: 1
              name: sample1
              value: 100
            - id: 2
              name: sample2
              value: 200
            - id: 3
              name: sample3
              value: 300
            '''.stripIndent()
        def cliParams = [samples: yamlFile.toString()]

        when:
        def samples = runScript(
            '''\
            params {
                samples: List<Record>
            }

            workflow { params.samples }
            ''',
            params: cliParams,
            configParams: [:]
        )

        then:
        samples instanceof List
        samples.size() == 3
        samples[0].id == 1
        samples[0].name == 'sample1'
        samples[0].value == 100
        samples[1].id == 2
        samples[2].id == 3

        cleanup:
        yamlFile?.delete()
    }

    def 'should load collection param from file specified in config'() {
        given:
        def jsonFile = Files.createTempFile('test', '.json')
        jsonFile.text = '[{"x": 1}, {"x": 2}]'
        def configParams = [items: jsonFile.toString()]

        when:
        def items = runScript(
            '''\
            params {
                items: List<Record>
            }

            workflow { params.items }
            ''',
            params: [:],
            configParams: configParams
        )

        then:
        items instanceof List
        items.size() == 2
        items[0].x == 1
        items[1].x == 2

        cleanup:
        jsonFile?.delete()
    }

    def 'should report error for unrecognized file format'() {
        given:
        def txtFile = Files.createTempFile('test', '.txt')
        txtFile.text = 'some text'
        def cliParams = [items: txtFile.toString()]

        when:
        runScript(
            '''\
            params {
                items: List<Record>
            }

            workflow { params.items }
            ''',
            params: cliParams,
            configParams: [:]
        )

        then:
        def e = thrown(ScriptRuntimeException)
        e.message.contains("Unrecognized file format 'txt'")
        e.message.contains("supplied for parameter `items` -- should be CSV, JSON, or YAML")

        cleanup:
        txtFile?.delete()
    }

    def 'should report error for invalid file content type'() {
        given:
        def jsonFile = Files.createTempFile('test', '.json')
        jsonFile.text = '{"not": "a list"}'
        def cliParams = [items: jsonFile.toString()]

        when:
        runScript(
            '''\
            params {
                items: List<Record>
            }

            workflow { params.items }
            ''',
            params: cliParams,
            configParams: [:]
        )

        then:
        def e = thrown(ScriptRuntimeException)
        e.message.contains('Parameter `items` with type List<Record> cannot be assigned to contents of')

        cleanup:
        jsonFile?.delete()
    }

    def 'should load record collection param from file'() {
        given:
        def inputFile = Files.createTempFile('test', '.json')
        inputFile.text = '''\
            [
              {"id": 1, "name": "sample1", "value": 100},
              {"id": 2, "name": "sample2", "value": 200},
              {"id": 3, "name": "sample3", "value": 300}
            ]
            '''.stripIndent()

        when:
        def samples = runScript(
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
            params: [samples: inputFile.toString()]
        )
        then:
        samples instanceof List
        samples.size() == 3
        samples[0] instanceof Record
        samples[0].id == 1
        samples[0].name == 'sample1'
        samples[0].value == 100
        samples[1].id == 2
        samples[2].id == 3
    }

    def 'should report error for missing record field'() {
        given:
        def inputFile = Files.createTempFile('test', '.json')
        inputFile.text = '''\
            [
              {"id": 1, "name": "sample1"}
            ]
            '''.stripIndent()

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
            params: [samples: inputFile.toString()]
        )

        then:
        def e = thrown(AbortOperationException)
        e.message == "Input record [id:1, name:sample1] is missing field 'value' required by record type 'Sample'"

        cleanup:
        inputFile?.delete()
    }

}
