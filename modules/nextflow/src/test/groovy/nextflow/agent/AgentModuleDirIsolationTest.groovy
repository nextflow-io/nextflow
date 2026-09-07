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
import nextflow.script.ScriptFile
import nextflow.script.ScriptRunner
import spock.lang.Timeout
import test.Dsl2Spec

/**
 * Regression test for the per-tool {@code moduleDir} bleed (Phase 3.x).
 *
 * <p>When an agent selects MULTIPLE processes from separate modules under
 * {@code nf:module_run}, each tool's {@code main.nf} must resolve {@code moduleDir} to ITS
 * OWN module directory. The agent used to compile a module path entry itself, as the
 * {@code mainScript} on the SHARED {@code session.binding}; {@code BaseScript.setup()}
 * overwrites {@code binding.moduleDir} per module, so the LAST-compiled module's dir won
 * (container/conda bleed across tools). That resolution form is gone — a module tool is now
 * {@code include}d like any other module, which is what gives each one its own isolated
 * {@link nextflow.script.ScriptBinding} — and this test pins that the property survived the
 * move: a bridge dispatching a tool must not disturb the module dir another tool resolves.
 *
 * <p>Three module tools live in SEPARATE temp dirs; each {@code exec}-only process
 * emits {@code moduleDir.toString()}. The test dispatches ALL THREE and asserts each tool
 * reports its OWN module dir and that they all DIFFER.
 *
 * <p>The {@code @Timeout} fails if the tool input queues are not poisoned on completion.
 */
@Timeout(90)
class AgentModuleDirIsolationTest extends Dsl2Spec {

    def cleanup() {
        AgentRunnerProvider.testRunner = null
    }

    def 'should resolve each external file-module tool to its own moduleDir'() {
        given:
        final root = Files.createTempDirectory('test')
        final work = root.resolve('work'); Files.createDirectories(work)

        // -- two distinct module dirs, each with a single scalar-`val`/`exec` process that
        //    emits its OWN moduleDir
        final modA = root.resolve('modA'); Files.createDirectories(modA)
        modA.resolve('main.nf').text = '''
            nextflow.enable.types = true

            process tool_a {
                input:
                x: String

                output:
                out: String

                exec:
                out = moduleDir.toString()
            }
            '''.stripIndent()
        final modADir = modA.toAbsolutePath().toString()

        final modB = root.resolve('modB'); Files.createDirectories(modB)
        modB.resolve('main.nf').text = '''
            nextflow.enable.types = true

            process tool_b {
                input:
                x: String

                output:
                out: String

                exec:
                out = moduleDir.toString()
            }
            '''.stripIndent()
        final modBDir = modB.toAbsolutePath().toString()

        final modC = root.resolve('modC'); Files.createDirectories(modC)
        modC.resolve('main.nf').text = '''
            nextflow.enable.types = true

            process tool_c {
                input:
                x: String

                output:
                out: String

                exec:
                out = moduleDir.toString()
            }
            '''.stripIndent()
        final modCDir = modC.toAbsolutePath().toString()

        and:
        final main = root.resolve('main.nf')
        main.text = """
            include { tool_a } from '${modA.resolve('main.nf').toAbsolutePath()}'
            include { tool_b } from '${modB.resolve('main.nf').toAbsolutePath()}'
            include { tool_c } from '${modC.resolve('main.nf').toAbsolutePath()}'

            agent a {
                model 'm'
                instruction 'i'
                tools 'nf:module_run:tool_a', 'nf:module_run:tool_b', 'nf:module_run:tool_c'

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
                a(channel.of('go')).view { it }
            }
            """.stripIndent()

        and:
        String resultA = null
        String resultB = null
        String resultC = null
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req ->
            // all three module tools are advertised
            assert (req.toolSpecs*.name as Set) == ['tool_a', 'tool_b', 'tool_c'] as Set
            // dispatch each tool: each runs its REAL process and emits its OWN moduleDir
            resultA = req.dispatch.call('tool_a', '{"x":"go"}')
            resultB = req.dispatch.call('tool_b', '{"x":"go"}')
            resultC = req.dispatch.call('tool_c', '{"x":"go"}')
            return resultA
        } as AgentRunner

        when:
        final runner = new ScriptRunner([process: [executor: 'local'], workDir: work.toString()])
        runner.setScript(new ScriptFile(main))
        runner.execute()

        then:
        resultA != null
        resultB != null
        resultC != null
        and:
        final outA = new JsonSlurper().parseText(resultA).out as String
        final outB = new JsonSlurper().parseText(resultB).out as String
        final outC = new JsonSlurper().parseText(resultC).out as String
        and:
        // each tool resolved its OWN module directory ...
        outA == modADir
        outB == modBDir
        outC == modCDir
        and:
        // ... and they are all distinct (no shared/last-compiled-module bleed)
        outA != outB
        outB != outC
        outA != outC
    }
}
