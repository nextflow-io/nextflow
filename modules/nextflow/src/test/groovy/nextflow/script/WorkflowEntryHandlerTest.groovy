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
            getDeclaredInputs() >> inputs
            getDeclaredInputTypes() >> [:]
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
            getDeclaredInputs() >> ['samples']
            getDeclaredInputTypes() >> [:]
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
        handler.getWorkflowArguments(workflowDef, [:])

        then:
        def e = thrown(ScriptRuntimeException)
        e.message.contains('requires input `samples`')
    }

    def 'should throw error when multiple workflows are defined'() {
        given:
        def workflow1 = Mock(WorkflowDef) {
            getName() >> 'FIRST'
            getDeclaredInputs() >> []
            getDeclaredInputTypes() >> [:]
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
            '''\
            nextflow.enable.types = true

            workflow GREET {
                take:
                name: String

                emit:
                greeting = "Hello, ${name}!"
            }
            ''',
            config: [params: [name: 'World']]
        )

        then:
        result != null
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
            '''\
            nextflow.enable.types = true

            workflow PROCESS_SAMPLES {
                take:
                samples: Channel<Record>

                emit:
                out = samples
            }
            ''',
            config: [params: [samples: csvFile.toString()]]
        )

        then:
        result != null

        cleanup:
        csvFile?.delete()
    }

    def 'should auto-run a named workflow with a JSON samplesheet input'() {
        given:
        def jsonFile = Files.createTempFile('samples', '.json')
        jsonFile.text = '[{"id":1,"name":"s1"},{"id":2,"name":"s2"}]'

        when:
        def result = runScript(
            '''\
            nextflow.enable.types = true

            workflow PROCESS_SAMPLES {
                take:
                samples: Channel<Record>

                emit:
                out = samples
            }
            ''',
            config: [params: [samples: jsonFile.toString()]]
        )

        then:
        result != null

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
            '''\
            nextflow.enable.types = true

            workflow PIPELINE {
                take:
                samples: Channel<Record>
                outdir: String

                emit:
                out = samples
            }
            ''',
            config: [params: [samples: csvFile.toString(), outdir: 'results']]
        )

        then:
        result != null

        cleanup:
        csvFile?.delete()
    }

    def 'should throw for a missing workflow input'() {
        when:
        runScript(
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

    def 'should prefer explicit entry workflow over named workflow'() {
        when:
        // An explicit (unnamed) entry workflow takes priority over WorkflowEntryHandler
        def result = runScript(
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

}
