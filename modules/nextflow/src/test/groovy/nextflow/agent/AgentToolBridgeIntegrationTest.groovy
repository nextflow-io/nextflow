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
import nextflow.exception.ScriptRuntimeException
import spock.lang.Timeout
import test.Dsl2Spec

import static test.ScriptHelper.runScript

/**
 * End-to-end test of the {@link ModuleToolBridge}: an agent declares an in-scope
 * process ({@code greet}) as a tool, and a mock runner invokes the dispatch
 * callback. This proves the headline mechanism — the LLM's tool call marshals
 * JSON args into channel values, the real {@code greet} process executes through
 * the standard dataflow/executor machinery, and its output is serialized back to
 * the caller as JSON — and that the run TERMINATES (the {@code @Timeout} fails if
 * the tool input queues are not poisoned on completion).
 */
@Timeout(60)
class AgentToolBridgeIntegrationTest extends Dsl2Spec {

    def cleanup() {
        AgentRunnerProvider.testRunner = null
    }

    def 'should run an in-scope process as an agent tool and terminate'() {
        given:
        AgentRunnerRequest captured = null
        String dispatchResult = null
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req ->
            captured = req
            // the bridge exposes a `greet` tool with a scalar `name:String` input
            assert req.toolSpecs.size() == 1
            assert req.toolSpecs[0].name == 'greet'
            assert req.toolSpecs[0].inputSchema.properties.name.type == 'string'
            // invoke the tool: this drives the REAL greet process through the executor
            dispatchResult = req.dispatch.call('greet', '{"name":"Ada"}')
            // the returned JSON proves the process actually ran and produced the value
            assert new JsonSlurper().parseText(dispatchResult) == [greeting: 'Hello Ada!']
            // the agent's final answer
            return dispatchResult
        } as AgentRunner

        when:
        def result = runScript('''
            nextflow.enable.types = true

            process greet {
                input:
                name: String

                output:
                greeting: String

                exec:
                greeting = "Hello ${name}!"
            }

            agent assistant {
                model 'm'
                instruction 'i'
                tools 'greet'

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
        // the workflow emits the runner's final answer (the dispatch result)
        new JsonSlurper().parseText(result.val) == [greeting: 'Hello Ada!']
        and:
        // the dispatch went through the bridge and returned the real process output
        captured != null
        new JsonSlurper().parseText(dispatchResult) == [greeting: 'Hello Ada!']
    }

    def 'should reject combining tools with a record (structured) output'() {
        given:
        // the runner must never be reached: the guard fires at agent run/build time
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req -> throw new IllegalStateException('runner should not be invoked') } as AgentRunner

        when:
        runScript('''
            nextflow.enable.types = true

            record Answer { greeting: String }

            process greet {
                input:
                name: String

                output:
                greeting: String

                exec:
                greeting = "Hello ${name}!"
            }

            agent assistant {
                model 'm'
                instruction 'i'
                tools 'greet'

                input:
                    request: String

                output:
                    answer: Answer

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
        def e = thrown(ScriptRuntimeException)
        e.message == 'Agent `assistant`: combining tools or skills with a record (structured) output is not yet supported - use a plain output type (e.g. String) when declaring tools or skills'
    }
}
