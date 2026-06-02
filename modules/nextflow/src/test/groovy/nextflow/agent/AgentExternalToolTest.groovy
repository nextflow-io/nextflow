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
import nextflow.exception.ScriptRuntimeException
import spock.lang.Timeout
import test.Dsl2Spec

import static test.ScriptHelper.runScript

/**
 * End-to-end test of an agent using an EXTERNAL module FILE as a tool (Phase 3.1).
 *
 * The {@code tools './mod.nf'} entry is compiled to a runnable {@link nextflow.script.ProcessDef}
 * on the owner's live session (pre-ignition) and pre-wired through the
 * {@link ModuleToolBridge}. A mock runner invokes the dispatch callback, proving the
 * external module's process actually executes through the standard dataflow/executor
 * machinery and its output is serialized back to the caller as JSON. The {@code @Timeout}
 * fails if the tool input queues are not poisoned on completion.
 */
@Timeout(60)
class AgentExternalToolTest extends Dsl2Spec {

    def cleanup() {
        AgentRunnerProvider.testRunner = null
    }

    private Path writeScripts(String mainScript, String moduleScript) {
        final dir = Files.createTempDirectory('test')
        dir.resolve('mod.nf').text = moduleScript.stripIndent()
        final main = dir.resolve('main.nf')
        main.text = mainScript.stripIndent()
        return main
    }

    def 'should compile an external module file as an agent tool and terminate'() {
        given:
        AgentRunnerRequest captured = null
        String dispatchResult = null
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req ->
            captured = req
            // the bridge exposes a `shout` tool with a scalar `text:String` input
            assert req.toolSpecs.size() == 1
            assert req.toolSpecs[0].name == 'shout'
            assert req.toolSpecs[0].inputSchema.properties.text.type == 'string'
            // invoke the tool: this drives the REAL external `shout` process through the executor
            dispatchResult = req.dispatch.call('shout', '{"text":"ada"}')
            // the returned JSON proves the external module's process actually ran
            assert new JsonSlurper().parseText(dispatchResult) == [result: 'ADA']
            // the agent's final answer
            return dispatchResult
        } as AgentRunner

        and:
        def main = writeScripts(
            '''
            nextflow.enable.types = true

            agent a {
                model 'm'
                instruction 'i'
                tools './mod.nf'

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
                a(channel.of('hi')).view { it }
            }
            ''',
            '''
            nextflow.enable.types = true

            process shout {
                input:
                text: String

                output:
                result: String

                exec:
                result = text.toUpperCase()
            }
            ''')

        when:
        def result = runScript(main)

        then:
        // the workflow emits the runner's final answer (the dispatch result)
        new JsonSlurper().parseText(result.val) == [result: 'ADA']
        and:
        // the dispatch went through the bridge and returned the real external process output
        captured != null
        new JsonSlurper().parseText(dispatchResult) == [result: 'ADA']
    }

    def 'should fail with a clear error when a registry module cannot be resolved (no registry access)'() {
        given:
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req -> throw new IllegalStateException('runner should not be invoked') } as AgentRunner

        and:
        // a registry ref that is NOT installed locally; with no registry access the resolution
        // fails and Phase 3.3 surfaces a clear, ref-naming error (no silent swallow)
        def main = writeScripts(
            '''
            nextflow.enable.types = true

            agent a {
                model 'm'
                instruction 'i'
                tools 'acme-bogus/does-not-exist'

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
                a(channel.of('hi')).view { it }
            }
            ''',
            'nextflow.enable.types = true')

        when:
        // point the registry at an unreachable URL so resolution fails fast (no real network)
        runScript([config: [registry: [url: 'http://127.0.0.1:1/api']]], main)

        then:
        def e = thrown(ScriptRuntimeException)
        e.message.contains('Unable to resolve agent tool module `acme-bogus/does-not-exist` from the registry')
        e.message.contains('Check registry access/credentials')
    }

    def 'should fail for a missing module file'() {
        given:
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req -> throw new IllegalStateException('runner should not be invoked') } as AgentRunner

        and:
        def main = writeScripts(
            '''
            nextflow.enable.types = true

            agent a {
                model 'm'
                instruction 'i'
                tools './missing.nf'

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
                a(channel.of('hi')).view { it }
            }
            ''',
            'nextflow.enable.types = true')

        when:
        runScript(main)

        then:
        def e = thrown(ScriptRuntimeException)
        e.message.contains('module file not found')
    }
}
