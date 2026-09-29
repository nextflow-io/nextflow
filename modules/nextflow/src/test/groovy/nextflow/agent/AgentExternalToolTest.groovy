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
 * The module is {@code include}d like any other, and the agent names its process under
 * {@code nf:module_run}; the resulting {@link nextflow.script.ProcessDef} is pre-wired through
 * the {@link ModuleToolBridge}. A mock runner invokes the dispatch callback, proving the
 * external module's process actually executes through the standard dataflow/executor
 * machinery and its output is serialized back to the caller as JSON. The {@code @Timeout}
 * fails if the tool input queues are not poisoned on completion.
 *
 * <p>The second test is the counterpart: the three reference shapes the directive used to
 * resolve by trial — a module path, a registry reference and a bare process name — are now
 * rejected by the grammar itself (G1), so the only way to a module tool is the {@code include}
 * the first test uses.
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

    def 'should run an included external module file as an agent tool and terminate'() {
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

            include { shout } from './mod.nf'

            agent a {
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

    def 'should reject the legacy #WHAT entry - the directive resolves no reference shape but a namespaced ref'() {
        given:
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req -> throw new IllegalStateException('runner should not be invoked') } as AgentRunner

        and:
        // NOTE the module file the harness writes DOES define `shout`, and `./missing.nf` does
        // NOT exist: neither fact matters any more. There is no resolution attempt to succeed or
        // fail, because the entry never gets past the grammar (G1) — which is precisely what
        // removing the fallthrough bought.
        def main = writeScripts(
            """
            nextflow.enable.types = true

            include { shout } from './mod.nf'

            agent a {
                model 'm'
                instruction 'i'
                tools '${ENTRY}'

                input:
                request: String

                output:
                answer: String

                prompt:
                \"\"\"
                \${request}
                \"\"\"
            }

            workflow {
                a(channel.of('hi')).view { it }
            }
            """,
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
        runScript(main)

        then:
        def e = thrown(ScriptRuntimeException)
        e.message.contains("Invalid tool reference `${ENTRY}`")
        and: 'the message points at the replacement rather than just refusing'
        e.message.contains('a tool reference must be namespaced as `family[:group]:name`')
        e.message.contains('`include` it and name its process')

        where:
        WHAT                | ENTRY
        'registry ref'      | 'acme-bogus/does-not-exist'
        'module path'       | './missing.nf'
        'bare process name' | 'shout'
    }
}
