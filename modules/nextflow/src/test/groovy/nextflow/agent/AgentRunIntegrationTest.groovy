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

import spock.lang.Timeout
import test.Dsl2Spec

import static test.ScriptHelper.runScript

/**
 * End-to-end test: a workflow invoking an {@code agent} runs as a dataflow
 * operator, rendering the prompt and delegating to the configured
 * {@link AgentRunner} (here a mock injected via the test seam).
 */
@Timeout(10)
class AgentRunIntegrationTest extends Dsl2Spec {

    def cleanup() {
        AgentRunnerProvider.testRunner = null
    }

    def 'should run a record-typed agent end-to-end against a mock runner'() {
        given:
        // capture the request the runner receives and echo a canned JSON answer
        AgentRunnerRequest captured = null
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req -> captured = req; '{"answer":"ok","confidence":0.9}' } as AgentRunner

        when:
        def result = runScript('''
            nextflow.enable.types = true

            record Question { text: String }
            record Answer { answer: String; confidence: Double }

            agent eval_agent {
                model 'openai/gpt-5-mini'
                instruction 'You are helpful.'
                tools()
                maxIterations 7

                input:
                    q: Question

                output:
                    a: Answer

                prompt:
                """
                Question: ${q.text}
                """
            }

            workflow {
                eval_agent(channel.of(record(text: 'analyze my reads')))
            }
            ''')

        then:
        // the emitted value is the bound output record
        def out = result.val
        out instanceof Map
        out.answer == 'ok'
        out.confidence == 0.9d
        and:
        captured.model == 'openai/gpt-5-mini'
        captured.instruction == 'You are helpful.'
        captured.maxIterations == 7
        captured.prompt.contains('Question: analyze my reads')
        and:
        // the output schema was derived from the Answer record type
        captured.outputSchema.properties.answer.type == 'string'
        captured.outputSchema.properties.confidence.type == 'number'
        and:
        // the input record was serialized to JSON
        captured.inputJson.contains('analyze my reads')
    }

    def 'should run a val-typed agent and emit the runner text verbatim'() {
        given:
        // a non-record (scalar) output opts out of structured output: the runner
        // returns plain text which is emitted verbatim (no JSON parse)
        AgentRunnerRequest captured = null
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req -> captured = req; 'hello world' } as AgentRunner

        when:
        def result = runScript('''
            nextflow.enable.types = true

            agent qa {
                model 'openai/gpt-5-mini'
                instruction 'You are helpful.'
                tools()

                input:
                    question: String

                output:
                    answer: String

                prompt:
                """
                Answer: ${question}
                """
            }

            workflow {
                qa(channel.of('what is FASTQ?'))
            }
            ''')

        then:
        // the emitted value is the runner text verbatim, NOT a parsed record
        result.val == 'hello world'
        and:
        // no structured-output schema was derived for a scalar output
        captured.outputSchema == null
        captured.prompt.contains('Answer: what is FASTQ?')
        and:
        // the scalar input was still serialized to JSON
        captured.inputJson.contains('what is FASTQ?')
    }
}
