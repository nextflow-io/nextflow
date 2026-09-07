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

import java.util.concurrent.CompletableFuture
import java.util.concurrent.Executors
import java.util.concurrent.TimeUnit

import groovy.json.JsonSlurper
import spock.lang.Timeout
import test.Dsl2Spec

import static test.ScriptHelper.runScript

/**
 * End-to-end tests (WITHOUT a real LLM) of MULTIPLE tools per agent and of
 * dispatch-level error handling.
 *
 * <p>The first test declares two in-scope processes — {@code shout} (uppercases)
 * and {@code whisper} (lowercases) — as tools. A mock runner asserts that BOTH
 * descriptors are advertised and then dispatches a call to EACH tool, proving the
 * two distinct processes ran (uppercase vs lowercase) and that the bridge routes
 * each call to the correct request-scoped process.
 *
 * <p>The second test exercises the dispatch-level error path: an unknown tool name
 * and an unparseable {@code argsJson} are each returned to the caller as a
 * {@code {"error": ...}} JSON string (so the LLM could recover) without the agent
 * loop crashing, and a subsequent valid call still runs the real process.
 *
 * <p>The {@code @Timeout} guards against the request gateway not being poisoned
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
                tools 'nf:module_run:shout', 'nf:module_run:whisper'

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

    def 'should return dispatch-level errors as JSON tool results without crashing the agent'() {
        given:
        AgentRunnerRequest captured = null
        String unknownToolResult = null
        String badJsonResult = null
        String okResult = null
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req ->
            captured = req
            // an unknown tool name -> a {"error": ...} naming the tool, NOT a thrown exception
            unknownToolResult = req.dispatch.call('nope', '{}')
            // unparseable args -> a {"error": ...} mentioning a parse failure
            badJsonResult = req.dispatch.call('shout', '{not json')
            // the bridge is still usable afterwards (the loop did not crash)
            okResult = req.dispatch.call('shout', '{"text":"Hi"}')
            return okResult
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

            agent assistant {
                model 'm'
                instruction 'i'
                tools 'nf:module_run:shout'

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
        captured != null
        and:
        // unknown-tool error is a well-formed JSON string naming the missing tool
        unknownToolResult instanceof String
        def unknown = new JsonSlurper().parseText(unknownToolResult)
        unknown.error != null
        unknown.error.contains('nope')
        and:
        // parse-failure error is a well-formed JSON string naming the tool + the parse problem
        badJsonResult instanceof String
        def bad = new JsonSlurper().parseText(badJsonResult)
        bad.error != null
        bad.error.contains('shout')
        bad.error.toLowerCase().contains('parse')
        and:
        // the agent loop did NOT crash: a subsequent valid call still runs the real process
        new JsonSlurper().parseText(okResult) == [loud: 'HI']
        new JsonSlurper().parseText(result.val) == [loud: 'HI']
    }

    def 'should correlate concurrent calls to the same tool when they complete out of order'() {
        given:
        final completionOrder = Collections.synchronizedList([])
        String slowResult
        String fastResult
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req ->
            final pool = Executors.newFixedThreadPool(2)
            try {
                final slow = CompletableFuture.supplyAsync({ ->
                    final value = req.dispatch.call('delayed', '{"text":"slow","delay":1200}')
                    completionOrder.add('slow')
                    return value
                }, pool)
                Thread.sleep(100)
                final fast = CompletableFuture.supplyAsync({ ->
                    final value = req.dispatch.call('delayed', '{"text":"fast","delay":0}')
                    completionOrder.add('fast')
                    return value
                }, pool)
                fastResult = fast.get(1, TimeUnit.SECONDS)
                slowResult = slow.get(5, TimeUnit.SECONDS)
                return fastResult
            }
            finally {
                pool.shutdownNow()
            }
        } as AgentRunner

        when:
        def result = runScript([config: [poolSize: 1]], '''
            nextflow.enable.types = true

            process delayed {
                input:
                text: String
                delay: Integer
                output:
                value: String
                exec:
                    Thread.sleep(delay)
                    value = text
            }

            agent assistant {
                model 'm'
                tools 'nf:module_run:delayed'
                input:
                request: String
                output:
                answer: String
                prompt: "${request}"
            }

            workflow {
                assistant(channel.of('run')).view { it }
            }
            ''')

        then:
        completionOrder == ['fast', 'slow']
        new JsonSlurper().parseText(fastResult) == [value: 'fast']
        new JsonSlurper().parseText(slowResult) == [value: 'slow']
        new JsonSlurper().parseText(result.val) == [value: 'fast']
    }

    def 'should keep a request-scoped invocation correlated across a process retry'() {
        given:
        String toolResult
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req ->
            toolResult = req.dispatch.call('flaky', '{"text":"recovered"}')
            return toolResult
        } as AgentRunner

        when:
        def result = runScript('''
            nextflow.enable.types = true

            process flaky {
                errorStrategy 'retry'
                maxRetries 1
                input:
                text: String
                output:
                value: String
                exec:
                    if( task.attempt == 1 )
                        throw new nextflow.exception.ProcessFailedException('fail first attempt')
                    value = text
            }

            agent assistant {
                model 'm'
                tools 'nf:module_run:flaky'
                input:
                request: String
                output:
                answer: String
                prompt: "${request}"
            }

            workflow {
                assistant(channel.of('run')).view { it }
            }
            ''')

        then:
        new JsonSlurper().parseText(toolResult) == [value: 'recovered']
        new JsonSlurper().parseText(result.val) == [value: 'recovered']
    }
}
