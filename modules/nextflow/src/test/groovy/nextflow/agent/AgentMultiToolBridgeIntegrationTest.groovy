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
import spock.lang.Timeout
import test.Dsl2Spec

import static test.ScriptHelper.runScript

/**
 * End-to-end test of MULTIPLE tools per agent, WITHOUT a real LLM.
 *
 * <p>The agent declares two in-scope processes — {@code shout} (uppercases) and
 * {@code whisper} (lowercases) — as tools. A mock runner asserts that BOTH
 * descriptors are advertised and then dispatches a call to EACH tool, proving the
 * two distinct processes ran (uppercase vs lowercase) and that the bridge routes
 * each call to the correct pre-wired process.
 *
 * <p>The {@code @Timeout} guards against the tool input queues not being poisoned
 * on completion (which would hang {@code session.await()}).
 */
@Timeout(60)
class AgentMultiToolBridgeIntegrationTest extends Dsl2Spec {

    def cleanup() {
        AgentRunnerProvider.testRunner = null
    }

    def 'should run multiple in-scope processes as agent tools and terminate'() {
        given:
        AgentRunnerRequest captured = null
        String shoutResult = null
        String whisperResult = null
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req ->
            captured = req
            // BOTH tool descriptors must be advertised to the LLM
            assert req.toolSpecs.size() == 2
            final names = req.toolSpecs*.name as Set
            assert names == ['shout', 'whisper'] as Set
            assert req.toolSpecs.find { it.name == 'shout' }.inputSchema.properties.text.type == 'string'
            assert req.toolSpecs.find { it.name == 'whisper' }.inputSchema.properties.text.type == 'string'
            // dispatch to EACH tool: this drives the two REAL processes through the executor
            shoutResult = req.dispatch.call('shout', '{"text":"Hi"}')
            whisperResult = req.dispatch.call('whisper', '{"text":"Hi"}')
            // distinct results prove the two distinct processes ran
            assert new JsonSlurper().parseText(shoutResult) == [loud: 'HI']
            assert new JsonSlurper().parseText(whisperResult) == [quiet: 'hi']
            return shoutResult
        } as AgentRunner

        when:
        def result = runScript('''
            nextflow.enable.types = true

            process shout {
                input:
                text: String

                output:
                loud: String

                exec:
                loud = text.toUpperCase()
            }

            process whisper {
                input:
                text: String

                output:
                quiet: String

                exec:
                quiet = text.toLowerCase()
            }

            agent assistant {
                model 'm'
                instruction 'i'
                tools 'shout', 'whisper'

                input:
                    request: String

                output:
                    answer: String

                prompt:
                """
                ${request}
                """
            }

            workflow {
                assistant(channel.of('hi')).view { it }
            }
            ''')

        then:
        // the workflow emits the runner's final answer (the shout dispatch result)
        new JsonSlurper().parseText(result.val) == [loud: 'HI']
        and:
        captured != null
        // both distinct processes ran and produced their distinct (upper/lower) outputs
        new JsonSlurper().parseText(shoutResult) == [loud: 'HI']
        new JsonSlurper().parseText(whisperResult) == [quiet: 'hi']
    }
}
