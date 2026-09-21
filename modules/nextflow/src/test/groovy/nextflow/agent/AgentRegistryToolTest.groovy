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

import groovy.json.JsonOutput
import groovy.json.JsonSlurper
import nextflow.module.ModuleChecksum
import nextflow.module.ModuleReference
import nextflow.module.ModuleStorage
import nextflow.script.ScriptFile
import nextflow.script.ScriptRunner
import spock.lang.Timeout
import test.Dsl2Spec

/**
 * End-to-end test of resolving a REGISTRY module reference (Phase 3.3) as an agent tool,
 * WITHOUT touching the network.
 *
 * The module is PRE-INSTALLED on disk in the layout {@link ModuleStorage} recognizes
 * ({@code <baseDir>/modules/<scope>/<name>/main.nf} + a sibling {@code meta.yml} manifest with
 * a {@code version}, plus a {@code .module-info} checksum so the install passes integrity). The
 * script {@code include}s the registry reference — {@code IncludeDef} builds a
 * {@link nextflow.module.ModuleResolver} rooted at {@code session.baseDir} and
 * {@code resolve(ref, null, autoInstall=true)} finds the local install WITHOUT any registry
 * call — and the agent then names the included process under {@code nf:module_run}. A registry
 * reference is no longer a {@code tools} entry of its own: the grammar has one way in, and it is
 * the same {@code include} every other module goes through. Since a sibling {@code meta.yml} is
 * present, schema/marshalling stays spec-driven (Phase 3.2).
 *
 * The module runs through the REAL local executor (it {@code cat}s the staged input into
 * {@code out.dat}), proving real resolution + staging + execution. The
 * {@code @Timeout} fails if the tool input queues are not poisoned on completion.
 */
@Timeout(90)
class AgentRegistryToolTest extends Dsl2Spec {

    def cleanup() {
        AgentRunnerProvider.testRunner = null
    }

    def 'should resolve a pre-installed registry module as an agent tool and run it (no network)'() {
        given:
        final dir = Files.createTempDirectory('test')
        final work = dir.resolve('work'); Files.createDirectories(work)
        final reads = dir.resolve('reads.txt'); reads.text = 'hello'
        final readsAbs = reads.toAbsolutePath().toString()

        and:
        // -- PRE-INSTALL the module locally so ModuleResolver.resolve returns WITHOUT network.
        //    Layout recognized by ModuleStorage.getInstalledModule:
        //      <baseDir>/modules/acme/echo/main.nf
        //      <baseDir>/modules/acme/echo/meta.yml   (manifest; MUST carry a `version`)
        //      <baseDir>/modules/acme/echo/.module-info (checksum -> integrity VALID)
        final reference = ModuleReference.parse('acme/echo')
        final storage = new ModuleStorage(dir)
        final moduleDir = storage.getModuleDir(reference)
        Files.createDirectories(moduleDir)

        // a CLASSIC DSL2 tuple module mimicking the nf-core shape; `cat` runs via the
        // LocalExecutor shell (no container)
        moduleDir.resolve('main.nf').text = '''
            process echo_tool {
                input:
                tuple val(meta), path(reads)

                output:
                tuple val(meta), path("out.dat"), emit: report

                script:
                """
                cat ${reads} > out.dat
                """
            }
            '''.stripIndent()

        // sibling meta.yml: the registry manifest + the tuple I/O spec (Phase 3.2)
        moduleDir.resolve('meta.yml').text = '''\
            name: acme/echo
            version: 1.0.0
            description: Copy the input reads to an output file
            input:
              - - name: meta
                  type: map
                  description: sample meta
                - name: reads
                  type: file
                  description: input reads
            output:
              - - name: meta
                  type: map
                  description: sample meta
                - name: outfile
                  type: file
                  description: the output
            '''.stripIndent()

        // mark the install as integrity-valid (checksum computed AFTER files are written)
        ModuleChecksum.save(moduleDir, ModuleChecksum.compute(moduleDir))

        and:
        final main = dir.resolve('main.nf')
        main.text = '''
            include { echo_tool } from 'acme/echo'

            agent a {
                model 'm'
                instruction 'i'
                tools 'nf:module_run:echo_tool'

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
                a(channel.of('go')).view { it }
            }
            '''.stripIndent()

        and:
        AgentRunnerRequest captured = null
        String dispatchResult = null
        Map<String,Object> outAssert = [:]
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req ->
            captured = req
            // the included registry module resolves to a single spec-driven tool whose
            // LLM-facing name is the PROCESS name — there is no reference-sanitizing step any
            // more, because the wire name comes from `nf:module_run:<PROCESS>` (§4)
            assert req.toolSpecs.size() == 1
            assert req.toolSpecs[0].name == 'echo_tool'
            assert req.toolSpecs[0].inputSchema.properties.meta.type == 'object'
            assert req.toolSpecs[0].inputSchema.properties.reads.type == 'string'
            // invoke the tool: drives the REAL echo_tool process through the executor,
            // staging `reads.txt` and producing `out.dat`
            dispatchResult = req.dispatch.call('echo_tool', JsonOutput.toJson([meta: [id: 's1'], reads: readsAbs]))
            final parsed = new JsonSlurper().parseText(dispatchResult) as Map
            // the result is keyed by the emit name `report`; its record carries `outfile`
            // as an ABSOLUTE PATH STRING and `meta` as the round-tripped object
            assert parsed.containsKey('report')
            final report = parsed.report as Map
            assert report.meta == [id: 's1']
            final outfile = report.outfile as String
            assert outfile.startsWith('/')
            final outPath = Path.of(outfile)
            assert Files.exists(outPath)
            assert outPath.text == 'hello'
            outAssert = [outfile: outfile, content: outPath.text]
            return dispatchResult
        } as AgentRunner

        when:
        // an unreachable registry URL keeps the test hermetic: the pre-installed module resolves
        // from disk, and the optional registry metadata fetch falls back to the sibling meta.yml
        final runner = new ScriptRunner([
            process: [executor: 'local'],
            workDir: work.toString(),
            registry: [url: 'http://localhost:0'] ])
        runner.setScript(new ScriptFile(main))
        runner.execute()

        then:
        // resolution + include + dispatch all succeeded and returned the real process output
        captured != null
        dispatchResult != null
        and:
        final parsed = new JsonSlurper().parseText(dispatchResult) as Map
        parsed.containsKey('report')
        (parsed.report as Map).meta == [id: 's1']
        ((parsed.report as Map).outfile as String).startsWith('/')
        and:
        // the file the real process produced exists and has the staged content
        outAssert.content == 'hello'
        Files.exists(Path.of(outAssert.outfile as String))
    }
}
