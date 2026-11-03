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

package nextflow.agent

import groovy.json.JsonSlurper
import nextflow.exception.AbortOperationException
import spock.lang.Timeout
import test.Dsl2Spec

import static test.ScriptHelper.runScript

/**
 * Tests for executing an agent directly via `nextflow module run`.
 */
@Timeout(30)
class AgentModuleRunTest extends Dsl2Spec {

    def cleanup() {
        AgentRunnerProvider.testRunner = null
    }

    def 'should execute a single agent with params as inputs'() {
        given:
        AgentRunnerRequest captured = null
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req ->
            captured = req; 'FASTQ is a text format'
        } as AgentRunner

        when:
        def result = runScript('''
            nextflow.enable.types = true

            agent qa {
                model 'openai/gpt-4o'
                input:
                question: String
                output:
                stdout()
                prompt: "Answer: ${question}"
            }
            ''', moduleRun: true, params: [question: 'What is FASTQ?'])

        then:
        captured.prompt == 'Answer: What is FASTQ?'
        result.val == 'FASTQ is a text format'
    }

    def 'should map params to a destructured record input'() {
        given:
        AgentRunnerRequest captured = null
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req ->
            captured = req; '{"score": 7}'
        } as AgentRunner

        when:
        def result = runScript('''
            nextflow.enable.types = true

            agent qa {
                model 'openai/gpt-4o'
                input:
                record(id: String, count: Integer)
                output:
                score: Integer
                prompt: "Score ${id} with ${count} reads"
            }
            ''', moduleRun: true, params: [id: 's1', count: '3'])

        then:
        new JsonSlurper().parseText(captured.inputJson) == [id: 's1', count: 3]
        result.val == 7
    }

    def 'should execute the agent rather than a process it uses as a tool'() {
        given:
        AgentRunnerRequest captured = null
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req ->
            captured = req; 'done'
        } as AgentRunner

        when:
        def result = runScript('''
            nextflow.enable.types = true

            process uppercase {
                input:
                text: String
                output:
                result: String
                exec:
                result = text.toUpperCase()
            }

            agent shouty {
                model 'openai/gpt-4o'
                tools 'nf:module_run'
                input:
                request: String
                output:
                stdout()
                prompt: "${request}"
            }
            ''', moduleRun: true, params: [request: 'shout hello'])

        then:
        captured.prompt == 'shout hello'
        captured.toolSpecs*.name == ['uppercase']
        result.val == 'done'
    }

    def 'should not execute an agent directly when it is not the only one, or without module run'() {
        given:
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req -> 'never' } as AgentRunner

        when:
        runScript(SCRIPT, moduleRun: MODULE_RUN, params: [q: 'hi'])

        then:
        def e = thrown(AbortOperationException)
        e.message.contains('No entry workflow specified')

        where:
        SCRIPT            | MODULE_RUN
        TWO_AGENTS        | true
        ONE_AGENT         | false
    }

    static final String ONE_AGENT = '''
        agent first {
            model 'openai/gpt-4o'
            input:
            q: String
            output:
            stdout()
            prompt: "${q}"
        }
        '''

    static final String TWO_AGENTS = ONE_AGENT + '''
        agent second {
            model 'openai/gpt-4o'
            input:
            q: String
            output:
            stdout()
            prompt: "${q}"
        }
        '''
}
