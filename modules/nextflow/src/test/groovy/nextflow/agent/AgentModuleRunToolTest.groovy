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
 * Integration test for the {@code module_run} capability: an agent declares
 * {@code tools 'module_run'} and all in-scope processes are enumerated and
 * surfaced to the LLM as a single {@code module_run} tool with a {@code module}
 * enum selector and a generic {@code args} object. A mock runner drives the
 * dispatch callback to verify end-to-end wiring and termination.
 */
@Timeout(60)
class AgentModuleRunToolTest extends Dsl2Spec {

    def cleanup() {
        AgentRunnerProvider.testRunner = null
    }

    def 'should expose in-scope processes as a single module_run tool'() {
        given:
        AgentRunnerRequest captured = null
        String dispatchResult = null
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req ->
            captured = req
            // must have exactly one tool: module_run
            assert req.toolSpecs.size() == 1
            assert req.toolSpecs[0].name == 'module_run'
            // the module enum must list the in-scope process name
            assert req.toolSpecs[0].inputSchema.properties.module.enum == ['greet']
            // args property must be present with additionalProperties=true
            assert req.toolSpecs[0].inputSchema.properties.args != null
            // invoke the tool: drives the REAL greet process through the executor
            dispatchResult = req.dispatch.call('module_run', '{"module":"greet","args":{"name":"Ada"}}')
            assert new JsonSlurper().parseText(dispatchResult) == [greeting: 'Hello Ada!']
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
                tools 'module_run'

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
        captured != null
        new JsonSlurper().parseText(dispatchResult) == [greeting: 'Hello Ada!']
    }

    def 'should return error for unknown module in module_run'() {
        given:
        String dispatchError = null
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req ->
            assert req.toolSpecs.size() == 1
            assert req.toolSpecs[0].name == 'module_run'
            // call with an unknown module name
            final result = req.dispatch.call('module_run', '{"module":"nope","args":{}}')
            dispatchError = result
            assert new JsonSlurper().parseText(result).error != null
            assert new JsonSlurper().parseText(result).error.contains('nope')
            return 'done'
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
                tools 'module_run'

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
        result.val == 'done'
        and:
        dispatchError != null
        new JsonSlurper().parseText(dispatchError).error.contains('nope')
    }

    def 'module_run and filesystem capabilities together produce two tool descriptors'() {
        given:
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req ->
            // must have two tools: module_run + filesystem
            assert req.toolSpecs.size() == 2
            assert req.toolSpecs*.name.contains('module_run')
            assert req.toolSpecs*.name.contains('filesystem')
            // invoke the module_run tool to verify dispatch still works
            final result = req.dispatch.call('module_run', '{"module":"greet","args":{"name":"World"}}')
            assert new JsonSlurper().parseText(result) == [greeting: 'Hello World!']
            return result
        } as AgentRunner

        when:
        runScript('''
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
                tools 'module_run', 'filesystem'

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
        noExceptionThrown()
    }
}
