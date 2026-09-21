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

import java.nio.file.Files
import java.nio.file.Path

import groovy.json.JsonSlurper
import spock.lang.TempDir
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

    @TempDir
    Path tempDir

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
                tools 'nf:module_run:greet'

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

    def 'should allow tools with a record (structured) output and bind the JSON'() {
        given:
        // the guard is gone (M5): a tool agent may also declare a structured output. The
        // plugin's final structuring turn is what returns schema JSON; here the stub runner
        // stands in for it and returns the JSON directly so the shared core bind runs.
        AgentRunnerRequest captured = null
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req ->
            captured = req
            return '{"greeting":"Hello Ada!"}'
        } as AgentRunner

        when:
        def result = runScript('''
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
                tools 'nf:module_run:greet'

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
        // the emitted value is a bound record with greeting == 'Hello Ada!'
        result.val.greeting == 'Hello Ada!'
        and:
        // the request carried a structured output schema (drives the plugin final-turn)
        captured != null
        captured.outputSchema != null
        captured.outputSchema.type == 'object'
        captured.outputSchema.properties.containsKey('greeting')
    }

    def 'a structured runner returning non-JSON aborts the run on the task path'() {
        given:
        // M-Tools: a tools agent now lowers to the task path, so a malformed structuring
        // answer fails the task body's JSON parse and aborts the run (the legacy path's clear
        // ScriptRuntimeException wrapper was removed together with runLegacy). The failure
        // names the agent/process and carries the JSON parse error in its cause chain.
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req -> return 'not json' } as AgentRunner

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
                tools 'nf:module_run:greet'

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
        // the malformed structuring answer fails the task body's JSON parse and aborts the run
        def e = thrown(Exception)
        hasCause(e, groovy.json.JsonException)
    }

    /** True when {@code type} appears anywhere in the throwable's cause chain. */
    private static boolean hasCause(Throwable t, Class type) {
        for( Throwable c = t; c != null; c = c.getCause() )
            if( type.isInstance(c) )
                return true
        return false
    }

    def 'should split a wrapper into N channels for a multi-output tools agent (PARTIAL)'() {
        given:
        // PARTIAL (M5): a tools agent with N>1 outputs requests a wrapper-object schema
        // and the wrapper is fanned out to one channel per output name.
        AgentRunnerRequest captured = null
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req ->
            captured = req
            return '{"plan":{"title":"Assembly"},"score":7}'
        } as AgentRunner

        when:
        def result = runScript('''
            nextflow.enable.types = true

            record Plan { title: String }

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
                tools 'nf:module_run:greet'

                input:
                request: String

                output:
                plan: Plan
                score: Integer

                prompt:
                """
                ${request}
                """
            }

            workflow {
                def r = assistant(channel.of('hi'))
                [ r.plan, r.score ]
            }
            ''')

        then:
        // the wrapper schema was requested: object root, both output names required
        captured != null
        captured.outputSchema.type == 'object'
        captured.outputSchema.required as Set == ['plan', 'score'] as Set
        captured.outputSchema.properties.keySet() as Set == ['plan', 'score'] as Set

        and:
        // the wrapper was split: each out.<name> channel received its slice, coerced to type
        def plan = result[0].val
        plan.title == 'Assembly'
        def score = result[1].val
        score == 7
    }

    def 'should allow skills with a record (structured) output and bind the JSON'() {
        given:
        // M5 decision 2: the guard is gone for skills too. A skills agent may declare a
        // structured output; the plugin's final structuring turn returns schema JSON (the
        // stub runner stands in for it here) and the SHARED core bind must parse+bind it.
        // This exercises the skills path through the exact same structured bind as the tools
        // case, at the core level (skills resolved on disk; runner stubbed).
        def skillDir = Files.createDirectories(tempDir.resolve('skills').resolve('greet'))
        Files.writeString(skillDir.resolve('SKILL.md'), "---\nname: greet\ndescription: greets\n---\ninstructions")
        def script = tempDir.resolve('main.nf')
        Files.writeString(script, '''
            nextflow.enable.types = true

            record Answer { greeting: String }

            agent assistant {
                model 'm'
                instruction 'i'
                skills 'greet'

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

        AgentRunnerRequest captured = null
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req ->
            captured = req
            return '{"greeting":"Hello Ada!"}'
        } as AgentRunner

        when:
        def result = runScript(script)

        then:
        // the emitted value is a bound record with greeting == 'Hello Ada!'
        result.val.greeting == 'Hello Ada!'
        and:
        // the skill resolved and the request carried it + the structured output schema
        captured != null
        captured.skills*.name == ['greet']
        captured.outputSchema != null
        captured.outputSchema.type == 'object'
        captured.outputSchema.properties.containsKey('greeting')
    }
}
