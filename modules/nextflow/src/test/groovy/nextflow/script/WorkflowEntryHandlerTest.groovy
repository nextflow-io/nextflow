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

import nextflow.Session
import nextflow.exception.AbortOperationException
import nextflow.exception.ScriptRuntimeException
import spock.lang.Timeout
import test.Dsl2Spec

import static test.ScriptHelper.*

/**
 * Tests for {@link WorkflowEntryHandler}.
 *
 * The per-type mapping of a param value to a declared input type is covered
 * by {@link ParamsHelperTest} -- the tests here only assert that a workflow
 * input routes through it.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Timeout(10)
class WorkflowEntryHandlerTest extends Dsl2Spec {

    // ── unit: loadFromFile ────────────────────────────────────────────────────

    private WorkflowEntryHandler makeHandler(List<String> inputs = [], List<String> workflowNames = ['HELLO']) {
        def workflowDef = Mock(WorkflowDef) {
            getName() >> workflowNames.first()
            getDeclaredInputs() >> inputs.collect { n -> new Param(n, null, false, null) }
        }
        def session = Mock(Session) { getParams() >> [:] }
        def script  = Mock(BaseScript) {
            isTypingEnabled() >> true
        }
        def meta    = Mock(ScriptMeta) {
            getLocalWorkflowNames() >> workflowNames
            getWorkflow(workflowNames.first()) >> workflowDef
        }
        return new WorkflowEntryHandler(script, session, meta)
    }

    def 'should load records from a samplesheet'() {
        given:
        def file = Files.createTempFile('test', ".${EXT}")
        file.text = TEXT

        when:
        def result = makeHandler().loadFromFile('samples', file.toAbsolutePath())

        then:
        result == EXPECTED

        cleanup:
        file?.delete()

        where:
        EXT    | TEXT                                             | EXPECTED
        // CSV has no types, so every value is a string
        'csv'  | 'id,name\n1,sample1\n2,sample2\n'                | [[id: '1', name: 'sample1'], [id: '2', name: 'sample2']]
        'json' | '[{"id":1,"name":"s1"},{"id":2,"name":"s2"}]'    | [[id: 1, name: 's1'], [id: 2, name: 's2']]
        'yml'  | '- id: 1\n  name: s1\n- id: 2\n  name: s2\n'     | [[id: 1, name: 's1'], [id: 2, name: 's2']]
    }

    def 'should throw for unrecognized samplesheet format'() {
        given:
        def txtFile = Files.createTempFile('test', '.txt')
        txtFile.text = 'some text'

        when:
        makeHandler().loadFromFile('items', txtFile.toAbsolutePath())

        then:
        def e = thrown(ScriptRuntimeException)
        e.message.contains("Unrecognized file format 'txt'")

        cleanup:
        txtFile?.delete()
    }

    def 'should throw for a JSON file whose top level is not a list'() {
        given:
        def jsonFile = Files.createTempFile('test', '.json')
        jsonFile.text = '{"key":"value"}'   // object, not array

        when:
        makeHandler().loadFromFile('samples', jsonFile.toAbsolutePath())

        then:
        def e = thrown(ScriptRuntimeException)
        e.message.contains('must contain a list of records')

        cleanup:
        jsonFile?.delete()
    }

    def 'should throw error when multiple workflows are defined'() {
        when:
        makeHandler([], ['FIRST', 'SECOND'])

        then:
        def e = thrown(IllegalStateException)
        e.message.contains('exactly one named workflow')
    }

    // ── integration tests ─────────────────────────────────────────────────────

    /**
     * Run a script as `nextflow module run`, with the given params supplied
     * either on the command line (raw strings) or in the config (structured
     * values) -- the two are resolved differently.
     */
    private static Object runWorkflow(String text, Map params, boolean fromConfig = false) {
        return fromConfig
            ? runScript(moduleRun: true, text, config: [params: params], configParams: params)
            : runScript(moduleRun: true, text, config: [params: params], params: params)
    }

    def 'should auto-run a named workflow with a scalar input'() {
        when:
        def result = runWorkflow(
            '''\
            nextflow.enable.types = true

            workflow GREET {
                take:
                name: String

                emit:
                greeting = "Hello, ${name}!"
            }
            ''',
            [name: 'World']
        )

        then:
        result.val == 'Hello, World!'
    }

    def 'should auto-run a named workflow with a CSV samplesheet input'() {
        given:
        def csvFile = Files.createTempFile('samples', '.csv')
        csvFile.text = '''\
            id,value
            1,alpha
            2,beta
            '''.stripIndent()

        when:
        def result = runWorkflow(
            '''\
            nextflow.enable.types = true

            workflow PROCESS_SAMPLES {
                take:
                samples: Channel<Record>

                emit:
                out = samples
            }
            ''',
            [samples: csvFile.toString()]
        )

        then:
        result.val == [id: '1', value: 'alpha']
        result.val == [id: '2', value: 'beta']

        cleanup:
        csvFile?.delete()
    }

    def 'should auto-run a named workflow with multiple inputs'() {
        given:
        def csvFile = Files.createTempFile('samples', '.csv')
        csvFile.text = '''\
            id,value
            1,alpha
            '''.stripIndent()

        when:
        def result = runWorkflow(
            '''\
            nextflow.enable.types = true

            workflow PIPELINE {
                take:
                samples: Channel<Record>
                outdir: String

                emit:
                out = samples.map { s -> "${outdir}/${s.id}" }
            }
            ''',
            [samples: csvFile.toString(), outdir: 'results']
        )

        then:
        result.val == 'results/1'

        cleanup:
        csvFile?.delete()
    }

    def 'should throw for a missing workflow input'() {
        when:
        runWorkflow(
            '''\
            nextflow.enable.types = true

            workflow GREET {
                take:
                name: String

                emit:
                greeting = "Hello!"
            }
            ''',
            [:]
        )

        then:
        def e = thrown(ScriptRuntimeException)
        e.message.contains('Parameter `--name` is required but no value was provided')
    }

    def 'should convert samplesheet records to the declared element type'() {
        given:
        def fastq = Files.createTempFile('reads', '.fastq')
        def csvFile = Files.createTempFile('samples', '.csv')
        csvFile.text = """\
            id,fastq_1,count,ratio
            s1,${fastq},7,1.5
            """.stripIndent()

        when:
        def result = runWorkflow(
            '''\
            nextflow.enable.types = true

            record Sample {
                id: String
                fastq_1: Path
                count: Integer
                ratio: Float
            }

            workflow RNASEQ {
                take:
                samples: Channel<Sample>

                emit:
                out = samples.map { s -> [s.getClass().simpleName, s.fastq_1.getClass().simpleName, s.count, s.ratio] }
            }
            ''',
            [samples: csvFile.toString()]
        )

        then:
        // each record is cast to the record type, and each field to its declared type
        result.val == ['RecordMap', 'UnixPath', 7, 1.5f]

        cleanup:
        csvFile?.delete()
        fastq?.delete()
    }

    def 'should reject a samplesheet that is missing a required record field'() {
        given:
        def csvFile = Files.createTempFile('samples', '.csv')
        csvFile.text = '''\
            id,wrong_column
            s1,oops
            '''.stripIndent()

        when:
        runWorkflow(
            '''\
            nextflow.enable.types = true

            record Sample {
                id: String
                fastq_1: Path
            }

            workflow RNASEQ {
                take:
                samples: Channel<Sample>

                emit:
                out = samples
            }
            ''',
            [samples: csvFile.toString()]
        )

        then:
        def e = thrown(Exception)
        e.message.contains("is missing field 'fastq_1' required by record type 'Sample'")

        cleanup:
        csvFile?.delete()
    }

    def 'should wrap a Value input in a channel'() {
        when:
        def result = runWorkflow(
            '''\
            nextflow.enable.types = true

            workflow CHECK {
                take:
                val: Value<Integer>

                emit:
                out = val.map { v -> [v.getClass().simpleName, v] }
            }
            ''',
            [val: '42']
        )

        then:
        // the param is resolved to the element type, then wrapped in a channel
        result.val == ['Integer', 42]
    }

    def 'should resolve a Path input and require it to exist'() {
        given:
        def file = Files.createTempFile('data', '.txt')

        when:
        def result = runWorkflow(
            '''\
            nextflow.enable.types = true

            workflow CHECK {
                take:
                data: Path

                emit:
                out = [data.getClass().simpleName, data.name]
            }
            ''',
            [data: file.toString()]
        )

        then:
        result.val == ['UnixPath', file.name]

        cleanup:
        file?.delete()
    }

    def 'should pass null for an omitted optional input'() {
        when:
        def result = runWorkflow(
            '''\
            nextflow.enable.types = true

            workflow OPT {
                take:
                required: String
                threads: Integer?

                emit:
                out = "${required}:${threads}"
            }
            ''',
            [required: 'hello']
        )

        then:
        result.val == 'hello:null'
    }

    def 'should reject a collection input given on the command line'() {
        when:
        runWorkflow(
            '''\
            nextflow.enable.types = true

            workflow COLL {
                take:
                ids: List<String>

                emit:
                out = ids
            }
            ''',
            [ids: 'a,b,c']
        )

        then:
        def e = thrown(ScriptRuntimeException)
        e.message.contains('cannot be assigned to a,b,c')
    }

    def 'should convert collection elements for an input from config'() {
        when:
        def result = runWorkflow(
            """\
            nextflow.enable.types = true

            workflow COLL {
                take:
                values: ${DECL}

                emit:
                out = [values.getClass().simpleName, values.toList().sort()]
            }
            """,
            [values: PARAM],
            true
        )

        then:
        result.val == [EXPECTED_CLASS, EXPECTED]

        where:
        DECL            | PARAM        | EXPECTED_CLASS  | EXPECTED
        'List<Integer>' | ['1', '12']  | 'ArrayList'     | [1, 12]
        'Set<Integer>'  | ['1', '2']   | 'LinkedHashSet' | [1, 2]
        'Bag<String>'   | ['a', 'b']   | 'HashBag'       | ['a', 'b']
    }

    def 'should reject a parameter that is not a workflow input'() {
        when:
        runWorkflow(
            '''\
            nextflow.enable.types = true

            workflow GREET {
                take:
                name: String

                emit:
                greeting = "Hello, ${name}!"
            }
            ''',
            [name: 'World', bogus: '1']
        )

        then:
        def e = thrown(ScriptRuntimeException)
        e.message.contains('Parameter `bogus` was specified on the command line but is not an input of workflow `GREET`')
    }

    def 'should prefer explicit entry workflow over named workflow'() {
        when:
        // An explicit (unnamed) entry workflow takes priority over WorkflowEntryHandler
        def result = runWorkflow(
            '''\
            workflow {
                "explicit entry"
            }

            workflow NAMED {
                take:
                x
                emit:
                out = x
            }
            ''',
            [x: 'ignored']
        )

        then:
        // the explicit entry workflow ran, and the named one was not invoked
        result == 'explicit entry'
    }

    def 'should reject a named workflow in a script that does not enable typing'() {
        when:
        runWorkflow(
            '''\
            workflow GREET {
                take:
                name
                emit:
                greeting = "Hello, ${name}!"
            }
            ''',
            [name: 'World']
        )

        then:
        def e = thrown(ScriptRuntimeException)
        e.message.contains('static typing is required')
    }

    def 'should reject a workflow input that is not typed'() {
        when:
        runWorkflow(
            '''\
            nextflow.enable.types = true

            workflow GREET {
                take:
                name

                emit:
                greeting = "Hello, ${name}!"
            }
            ''',
            [name: 'World']
        )

        then:
        def e = thrown(ScriptRuntimeException)
        e.message.contains('the following inputs are not typed: name')
    }

    def 'should report a samplesheet file that does not exist'() {
        when:
        runWorkflow(
            '''\
            nextflow.enable.types = true

            workflow RNASEQ {
                take:
                samples: Channel<Record>

                emit:
                out = samples
            }
            ''',
            [samples: '/some/missing/samples.csv']
        )

        then:
        // the file is reported by name, rather than as a raw NoSuchFileException
        def e = thrown(Exception)
        e.message.contains("Input file '/some/missing/samples.csv' does not exist")
    }

    def 'should not execute a named workflow directly without module run'() {
        when:
        runScript(
            '''\
            nextflow.enable.types = true

            workflow GREET {
                take:
                name: String

                emit:
                greeting = "Hello, ${name}!"
            }
            ''',
            config: [params: [name: 'World']],
            params: [name: 'World']
        )

        then:
        def e = thrown(AbortOperationException)
        e.message.contains('No entry workflow specified')
    }

    def 'should report a required input with no value'() {
        when:
        // e.g. a blank entry in a params file -- a required input must not be
        // silently satisfied by null
        runWorkflow(
            '''\
            nextflow.enable.types = true

            workflow GREET {
                take:
                name: String

                emit:
                greeting = "Hello, ${name}!"
            }
            ''',
            [name: null]
        )

        then:
        def e = thrown(ScriptRuntimeException)
        e.message.contains('Parameter `--name` is required but no value was provided')
    }

    def 'should pass a null param value to a nullable input'() {
        when:
        def result = runWorkflow(
            '''\
            nextflow.enable.types = true

            workflow GREET {
                take:
                name: String?

                emit:
                greeting = "Hello, ${name ?: 'World'}!"
            }
            ''',
            [name: null]
        )

        then:
        result.val == 'Hello, World!'
    }

    def 'should reject a channel element type that is not a record'() {
        given:
        def file = Files.createTempFile('samples', '.csv')
        file.text = 'id\ns1\ns2\n'

        when:
        // a samplesheet record would otherwise be stringified, e.g. '[id:s1]'
        runWorkflow(
            '''\
            nextflow.enable.types = true

            workflow GREET {
                take:
                samples: Channel<String>

                emit:
                out = samples
            }
            ''',
            [samples: file.toString()]
        )

        then:
        def e = thrown(ScriptRuntimeException)
        e.message.contains('cannot be loaded from a samplesheet')

        cleanup:
        file?.delete()
    }

    def 'should report an invalid samplesheet record with the param and file'() {
        given:
        def file = Files.createTempFile('samples', '.csv')
        file.text = 'id,count\ns1,abc\n'

        when:
        runWorkflow(
            '''\
            nextflow.enable.types = true

            record Sample {
                id: String
                count: Integer
            }

            workflow GREET {
                take:
                samples: Channel<Sample>

                emit:
                out = samples
            }
            ''',
            [samples: file.toString()]
        )

        then: 'the param and file are named, rather than a bare NumberFormatException'
        def e = thrown(ScriptRuntimeException)
        e.message.contains('Invalid record in samplesheet')
        e.message.contains('workflow input `samples`')

        cleanup:
        file?.delete()
    }

    def 'should report a collection element that cannot be converted'() {
        when:
        runWorkflow(
            '''\
            nextflow.enable.types = true

            workflow GREET {
                take:
                ids: List<Integer>

                emit:
                out = ids
            }
            ''',
            [ids: ['a', 'b']],
            true
        )

        then: 'the param and declared type are named, rather than a bare NumberFormatException'
        def e = thrown(ScriptRuntimeException)
        e.message.contains('Parameter `ids` with type List<Integer> cannot be assigned to')
    }

    def 'should load a blank samplesheet cell as a missing value'() {
        given:
        def fastq = Files.createTempFile('reads', '.fq')
        def file = Files.createTempFile('samples', '.csv')
        file.text = "id,fastq_1,fastq_2\ns1,${fastq},\n"

        when: 'an optional field is left blank rather than omitted'
        def result = runWorkflow(
            '''\
            nextflow.enable.types = true

            record Sample {
                id: String
                fastq_1: Path
                fastq_2: Path?
            }

            workflow GREET {
                take:
                samples: Channel<Sample>

                emit:
                out = samples.map { s -> "${s.id}:${s.fastq_2}" }
            }
            ''',
            [samples: file.toString()]
        )

        then:
        result.val == 's1:null'

        cleanup:
        file?.delete()
        fastq?.delete()
    }

    def 'should report a required record field left blank in a samplesheet'() {
        given:
        def file = Files.createTempFile('samples', '.csv')
        file.text = 'id,count\ns1,\n'

        when:
        runWorkflow(
            '''\
            nextflow.enable.types = true

            record Sample {
                id: String
                count: Integer
            }

            workflow GREET {
                take:
                samples: Channel<Sample>

                emit:
                out = samples
            }
            ''',
            [samples: file.toString()]
        )

        then: 'the missing field is named, rather than an empty-string conversion failure'
        def e = thrown(Exception)
        e.message.contains("missing field 'count'")

        cleanup:
        file?.delete()
    }

}
