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
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Timeout(10)
class WorkflowEntryHandlerTest extends Dsl2Spec {

    // ── unit: loadFromFile ────────────────────────────────────────────────────

    private WorkflowEntryHandler makeHandler(List<String> inputs = []) {
        def workflowDef = Mock(WorkflowDef) {
            getName() >> 'HELLO'
            getDeclaredInputs() >> inputs.collect { n -> new Param(n, null, false, null) }
        }
        def session = Mock(Session) { getParams() >> [:] }
        def script  = Mock(BaseScript) {
            isTypingEnabled() >> true
        }
        def meta    = Mock(ScriptMeta) {
            getLocalWorkflowNames() >> ['HELLO']
            getWorkflow('HELLO') >> workflowDef
        }
        return new WorkflowEntryHandler(script, session, meta)
    }

    def 'should load records from a CSV file'() {
        given:
        def csvFile = Files.createTempFile('test', '.csv')
        csvFile.text = '''\
            id,name
            1,sample1
            2,sample2
            '''.stripIndent()

        when:
        def result = makeHandler().loadFromFile('samples', csvFile.toAbsolutePath())

        then:
        result instanceof List
        result.size() == 2
        result[0].id == '1'
        result[0].name == 'sample1'

        cleanup:
        csvFile?.delete()
    }

    def 'should load records from a JSON file'() {
        given:
        def jsonFile = Files.createTempFile('test', '.json')
        jsonFile.text = '[{"id":1,"name":"s1"},{"id":2,"name":"s2"}]'

        when:
        def result = makeHandler().loadFromFile('samples', jsonFile.toAbsolutePath())

        then:
        result instanceof List
        result.size() == 2
        result[0].id == 1
        result[1].name == 's2'

        cleanup:
        jsonFile?.delete()
    }

    def 'should load records from a YAML file'() {
        given:
        def yamlFile = Files.createTempFile('test', '.yml')
        yamlFile.text = '''\
            - id: 1
              name: s1
            - id: 2
              name: s2
            '''.stripIndent()

        when:
        def result = makeHandler().loadFromFile('samples', yamlFile.toAbsolutePath())

        then:
        result instanceof List
        result.size() == 2
        result[0].id == 1
        result[1].name == 's2'

        cleanup:
        yamlFile?.delete()
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

    // ── unit: getWorkflowArguments / error cases ──────────────────────────────

    def 'should throw for a missing required workflow input'() {
        given:
        def workflowDef = Mock(WorkflowDef) {
            getName() >> 'HELLO'
            getDeclaredInputs() >> [new Param('samples', String, false, null)]
        }
        def session = Mock(Session)
        def script  = Mock(BaseScript) {
            isTypingEnabled() >> true
        }
        def meta    = Mock(ScriptMeta) {
            getLocalWorkflowNames() >> ['HELLO']
            getWorkflow('HELLO') >> workflowDef
        }
        def handler = new WorkflowEntryHandler(script, session, meta)

        when:
        handler.getWorkflowArguments(workflowDef)

        then:
        def e = thrown(ScriptRuntimeException)
        e.message.contains('requires input `samples`')
    }

    def 'should throw error when multiple workflows are defined'() {
        given:
        def workflow1 = Mock(WorkflowDef) {
            getName() >> 'FIRST'
            getDeclaredInputs() >> []
        }
        def session = Mock(Session) { getParams() >> [:] }
        def script  = Mock(BaseScript) {
            isTypingEnabled() >> true
        }
        def meta    = Mock(ScriptMeta) {
            getLocalWorkflowNames() >> ['FIRST', 'SECOND']
            getWorkflow('FIRST') >> workflow1
        }

        when:
        def handler = new WorkflowEntryHandler(script, session, meta)

        then:
        def e = thrown(IllegalStateException)
        e.message.contains('exactly one named workflow')
    }

    // ── integration tests ─────────────────────────────────────────────────────

    def 'should auto-run a named workflow with a scalar input'() {
        when:
        def result = runScript(
            moduleRun: true,
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
        def result = runScript(
            moduleRun: true,
            '''\
            nextflow.enable.types = true

            workflow PROCESS_SAMPLES {
                take:
                samples: Channel<Record>

                emit:
                out = samples
            }
            ''',
            config: [params: [samples: csvFile.toString()]],
            params: [samples: csvFile.toString()]
        )

        then:
        result.val == [id: '1', value: 'alpha']
        result.val == [id: '2', value: 'beta']

        cleanup:
        csvFile?.delete()
    }

    def 'should auto-run a named workflow with a JSON samplesheet input'() {
        given:
        def jsonFile = Files.createTempFile('samples', '.json')
        jsonFile.text = '[{"id":1,"name":"s1"},{"id":2,"name":"s2"}]'

        when:
        def result = runScript(
            moduleRun: true,
            '''\
            nextflow.enable.types = true

            workflow PROCESS_SAMPLES {
                take:
                samples: Channel<Record>

                emit:
                out = samples
            }
            ''',
            config: [params: [samples: jsonFile.toString()]],
            params: [samples: jsonFile.toString()]
        )

        then:
        result.val == [id: 1, name: 's1']
        result.val == [id: 2, name: 's2']

        cleanup:
        jsonFile?.delete()
    }

    def 'should auto-run a named workflow with multiple inputs'() {
        given:
        def csvFile = Files.createTempFile('samples', '.csv')
        csvFile.text = '''\
            id,value
            1,alpha
            '''.stripIndent()

        when:
        def result = runScript(
            moduleRun: true,
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
            config: [params: [samples: csvFile.toString(), outdir: 'results']],
            params: [samples: csvFile.toString(), outdir: 'results']
        )

        then:
        result.val == 'results/1'

        cleanup:
        csvFile?.delete()
    }

    def 'should throw for a missing workflow input'() {
        when:
        runScript(
            moduleRun: true,
            '''\
            nextflow.enable.types = true

            workflow GREET {
                take:
                name: String

                emit:
                greeting = "Hello!"
            }
            ''',
            params: [:]
        )

        then:
        def e = thrown(ScriptRuntimeException)
        e.message.contains('requires input `name`')
    }

    def 'should convert samplesheet records to the declared element type'() {
        given:
        def fastq = Files.createTempFile('reads', '.fastq')
        def csvFile = Files.createTempFile('samples', '.csv')
        csvFile.text = """\
            id,fastq_1
            s1,${fastq}
            """.stripIndent()

        when:
        def result = runScript(
            moduleRun: true,
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
                out = samples.map { s -> [s.getClass().simpleName, s.fastq_1.getClass().simpleName] }
            }
            ''',
            config: [params: [samples: csvFile.toString()]],
            params: [samples: csvFile.toString()]
        )

        then:
        // each record is cast to the record type and `fastq_1` to a Path
        result.val == ['RecordMap', 'UnixPath']

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
        runScript(
            moduleRun: true,
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
            config: [params: [samples: csvFile.toString()]],
            params: [samples: csvFile.toString()]
        )

        then:
        def e = thrown(Exception)
        e.message.contains("is missing field 'fastq_1' required by record type 'Sample'")

        cleanup:
        csvFile?.delete()
    }

    def 'should coerce scalar inputs to the declared type'() {
        when:
        def params = [dur: '2h', mem: '4.GB', ver: '1.2.3', num: '5.5', val: '42']
        def result = runScript(
            moduleRun: true,
            '''\
            nextflow.enable.types = true

            workflow CHECK {
                take:
                dur: Duration
                mem: MemoryUnit
                ver: VersionNumber
                num: Float
                val: Value<Integer>

                emit:
                out = val.map { v -> [
                    dur.getClass().simpleName,
                    mem.getClass().simpleName,
                    ver.getClass().simpleName,
                    num.getClass().simpleName,
                    v.getClass().simpleName
                ] }
            }
            ''',
            config: [params: params],
            params: params
        )

        then:
        result.val == ['Duration', 'MemoryUnit', 'VersionNumber', 'Float', 'Integer']
    }

    def 'should pass null for an omitted optional input'() {
        when:
        def result = runScript(
            moduleRun: true,
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
            config: [params: [required: 'hello']],
            params: [required: 'hello']
        )

        then:
        result.val == 'hello:null'
    }

    def 'should reject a collection input given on the command line'() {
        when:
        runScript(
            moduleRun: true,
            '''\
            nextflow.enable.types = true

            workflow COLL {
                take:
                ids: List<String>

                emit:
                out = ids
            }
            ''',
            config: [params: [ids: 'a,b,c']],
            params: [ids: 'a,b,c']
        )

        then:
        def e = thrown(ScriptRuntimeException)
        e.message.contains('cannot be assigned to a,b,c')
    }

    def 'should convert collection elements for an input from config'() {
        when:
        def result = runScript(
            moduleRun: true,
            '''\
            nextflow.enable.types = true

            workflow COLL {
                take:
                sizes: List<Integer>

                emit:
                out = sizes
            }
            ''',
            config: [params: [sizes: ['1', '12']]],
            configParams: [sizes: ['1', '12']]
        )

        then:
        result.val == [1, 12]
    }

    def 'should reject a parameter that is not a workflow input'() {
        when:
        runScript(
            moduleRun: true,
            '''\
            nextflow.enable.types = true

            workflow GREET {
                take:
                name: String

                emit:
                greeting = "Hello, ${name}!"
            }
            ''',
            config: [params: [name: 'World', bogus: '1']],
            params: [name: 'World', bogus: '1']
        )

        then:
        def e = thrown(ScriptRuntimeException)
        e.message.contains('Parameter `bogus` was specified on the command line but is not an input of workflow `GREET`')
    }

    def 'should prefer explicit entry workflow over named workflow'() {
        when:
        // An explicit (unnamed) entry workflow takes priority over WorkflowEntryHandler
        def result = runScript(
            moduleRun: true,
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
            params: [x: 'ignored']
        )

        then:
        // The explicit entry workflow ran
        result != null
    }

    def 'should coerce an integral command line value to a Float input'() {
        when:
        def result = runScript(
            moduleRun: true,
            '''\
            nextflow.enable.types = true

            workflow CHECK {
                take:
                num: Float

                emit:
                out = [num.getClass().simpleName, num]
            }
            ''',
            config: [params: [num: '5']],
            params: [num: '5']
        )

        then:
        result.val == ['Float', 5.0f]
    }

    def 'should resolve a Path input and require it to exist'() {
        given:
        def file = Files.createTempFile('data', '.txt')

        when:
        def result = runScript(
            moduleRun: true,
            '''\
            nextflow.enable.types = true

            workflow CHECK {
                take:
                data: Path

                emit:
                out = [data.getClass().simpleName, data.name]
            }
            ''',
            config: [params: [data: file.toString()]],
            params: [data: file.toString()]
        )

        then:
        result.val == ['UnixPath', file.name]

        cleanup:
        file?.delete()
    }

    def 'should reject a Path input that does not exist'() {
        when:
        runScript(
            moduleRun: true,
            '''\
            nextflow.enable.types = true

            workflow CHECK {
                take:
                data: Path

                emit:
                out = data
            }
            ''',
            config: [params: [data: '/some/missing/file.txt']],
            params: [data: '/some/missing/file.txt']
        )

        then:
        def e = thrown(Exception)
        e.message.contains('/some/missing/file.txt')
    }

    def 'should convert collection elements for Set and Bag inputs from config'() {
        when:
        def params = [ids: ['1', '2'], tags: ['a', 'b']]
        def result = runScript(
            moduleRun: true,
            '''\
            nextflow.enable.types = true

            workflow CHECK {
                take:
                ids: Set<Integer>
                tags: Bag<String>

                emit:
                out = [ids.getClass().simpleName, ids.toList(), tags.getClass().simpleName, tags.size()]
            }
            ''',
            config: [params: params],
            configParams: params
        )

        then:
        result.val == ['LinkedHashSet', [1, 2], 'HashBag', 2]
    }

    def 'should convert numeric record fields from a samplesheet'() {
        given:
        def csvFile = Files.createTempFile('samples', '.csv')
        csvFile.text = '''\
            id,count,ratio
            s1,7,1.5
            '''.stripIndent()

        when:
        def result = runScript(
            moduleRun: true,
            '''\
            nextflow.enable.types = true

            record Sample {
                id: String
                count: Integer
                ratio: Float
            }

            workflow RNASEQ {
                take:
                samples: Channel<Sample>

                emit:
                out = samples.map { s -> [s.count, s.ratio] }
            }
            ''',
            config: [params: [samples: csvFile.toString()]],
            params: [samples: csvFile.toString()]
        )

        then:
        result.val == [7, 1.5f]

        cleanup:
        csvFile?.delete()
    }

    def 'should reject a named workflow in a script that does not enable typing'() {
        when:
        runScript(
            moduleRun: true,
            '''\
            workflow GREET {
                take:
                name
                emit:
                greeting = "Hello, ${name}!"
            }
            ''',
            params: [name: 'World']
        )

        then:
        def e = thrown(ScriptRuntimeException)
        e.message.contains('does not enable static typing')
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

}
