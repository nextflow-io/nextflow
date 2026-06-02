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

    def 'should run an agent end-to-end against a mock runner'() {
        given:
        // capture the request the runner receives and echo a canned answer
        AgentRunnerRequest captured = null
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req -> captured = req; "PLAN: do fastqc" } as AgentRunner

        when:
        def result = runScript('''
            nextflow.enable.types = true

            agent eval_agent {
                model 'openai/gpt-5-mini'
                instruction 'You are helpful.'
                tools()
                maxIterations 7

                input:
                    question: String

                output:
                    plan: String

                prompt:
                """
                Question: ${question}
                """
            }

            workflow {
                eval_agent(channel.of('analyze my reads'))
            }
            ''')

        then:
        result.val == 'PLAN: do fastqc'
        and:
        captured.model == 'openai/gpt-5-mini'
        captured.instruction == 'You are helpful.'
        captured.maxIterations == 7
        captured.prompt.contains('Question: analyze my reads')
    }
}
