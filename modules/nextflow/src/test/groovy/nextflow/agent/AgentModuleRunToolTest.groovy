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
import nextflow.script.ScriptFile
import nextflow.script.ScriptRunner
import spock.lang.Timeout
import test.Dsl2Spec

import static test.ScriptHelper.runScript

/**
 * Integration test for the {@code module_run} capability: an agent declares
 * {@code tools 'module_run'} and every in-scope / included process is exposed to the LLM
 * as its OWN tool whose function {@code parameters} schema IS that module's flattened input
 * schema (required fields, {@code additionalProperties:false}, the nf-core {@code meta.id}
 * convention). This per-module enforcement is what lets OpenAI function-calling validate the
 * field names against the module schema — a single aggregate {@code module_run} tool could only
 * use a generic {@code args:{additionalProperties:true}} object, which OpenAI cannot enforce.
 *
 * A mock runner drives the dispatch callback to verify end-to-end wiring and termination.
 */
@Timeout(60)
class AgentModuleRunToolTest extends Dsl2Spec {

    def cleanup() {
        AgentRunnerProvider.testRunner = null
    }

    def 'should expose an in-scope process as its own per-module tool (not module_run)'() {
        given:
        AgentRunnerRequest captured = null
        String dispatchResult = null
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req ->
            captured = req
            // exactly one tool, NAMED after the module — NOT an aggregate 'module_run'
            assert req.toolSpecs.size() == 1
            assert req.toolSpecs[0].name == 'greet'
            assert req.toolSpecs*.name.every { it != 'module_run' }
            // the tool's parameters schema IS the module's flattened input schema
            assert req.toolSpecs[0].inputSchema.properties.name != null
            // dispatch by the MODULE tool name (not 'module_run'): drives the REAL greet process
            dispatchResult = req.dispatch.call('greet', '{"name":"Ada"}')
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

    def 'should return error for an unknown tool name'() {
        given:
        String dispatchError = null
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req ->
            assert req.toolSpecs.size() == 1
            assert req.toolSpecs[0].name == 'greet'
            // call with an unknown tool name
            final result = req.dispatch.call('nope', '{}')
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

    @Timeout(120)
    def 'should expose an included classic module (meta.yml-driven) as its own tool with an ENFORCED flattened schema, run it, and whitelist its output dir'() {
        given:
        final dir = Files.createTempDirectory('test')
        final work = dir.resolve('work'); Files.createDirectories(work)
        final moduleDir = dir.resolve('module'); Files.createDirectories(moduleDir)
        final reads = dir.resolve('reads.txt'); reads.text = 'hello'
        final readsAbs = reads.toAbsolutePath().toString()

        // a CLASSIC untyped module mimicking the nf-core tuple shape; the input file param is
        // named `fastq` (the convention the model must NOT rename/omit). The script `cat`s the
        // staged input into out.dat (a NON text-like extension, so the output stays an opaque
        // absolute path handle — exercising the output-dir whitelist relocation).
        moduleDir.resolve('main.nf').text = '''\
            process FOO {
                input:
                tuple val(meta), path(fastq)

                output:
                tuple val(meta), path('out.dat'), emit: report

                script:
                """
                cat ${fastq} > out.dat
                """
            }
            '''.stripIndent()

        // sibling meta.yml describing the tuple I/O (nf-core style); `fastq` is a required file
        moduleDir.resolve('meta.yml').text = '''\
            name: FOO
            description: Classic nf-core style module
            input:
              - - name: meta
                  type: map
                  description: sample meta
                - name: fastq
                  type: file
                  description: input reads
            output:
              - - name: meta
                  type: map
                  description: sample meta
                - name: outfile
                  type: file
                  description: output file
            '''.stripIndent()

        final main = dir.resolve('main.nf')
        main.text = """\
            include { FOO } from './module/main.nf'

            agent assistant {
                model 'm'
                instruction 'i'
                tools 'module_run', 'filesystem'

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
                assistant(channel.of('hi')).view { it }
            }
            """.stripIndent()

        and:
        AgentRunnerRequest captured = null
        String dispatchResult = null
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req ->
            captured = req
            // exactly two tools: the per-module tool FOO + filesystem — NOT an aggregate 'module_run'
            assert req.toolSpecs.size() == 2
            assert req.toolSpecs*.name.contains('FOO')
            assert req.toolSpecs*.name.contains('filesystem')
            assert req.toolSpecs*.name.every { it != 'module_run' }
            final fooSpec = req.toolSpecs.find { it.name == 'FOO' }
            assert req.toolSpecs*.name.every { it != 'module_run' }

            // the tool's parameters schema IS the module's flattened input schema, ENFORCED:
            //     `fastq` is a property, `required` includes the file input, and
            //     additionalProperties is false; `meta` is a nested object with the meta.id convention
            final schema = fooSpec.inputSchema
            final props = schema.properties as Map
            assert props.containsKey('fastq')
            assert props.containsKey('meta')
            // file input flattened to a string property
            assert (props.fastq as Map).type == 'string'
            // the file input is required (the enforcement the LLM cannot omit)
            assert (schema.required as List).contains('fastq')
            // OpenAI function-calling forbids extra/renamed fields
            assert schema.additionalProperties == false
            // nf-core meta.id convention: meta is a nested object carrying an `id`
            final meta = props.meta as Map
            assert meta.type == 'object'
            assert (meta.properties as Map).containsKey('id')

            // (c) dispatch by the module's OWN tool name runs the REAL process and returns output
            dispatchResult = req.dispatch.call('FOO', JsonOutput.toJson([meta: [id: 's1'], fastq: readsAbs]))
            final parsed = new JsonSlurper().parseText(dispatchResult) as Map
            assert parsed.containsKey('report')
            final report = parsed.report as Map
            assert report.meta == [id: 's1']
            final outfile = report.outfile as String
            assert outfile.startsWith('/')

            // (d) after the run, the dispatch context's readable dirs include the module output's
            //     parent dir (whitelist relocation): a subsequent filesystem read would be allowed
            final ctx = getDispatchContext()
            assert ctx != null
            final outParent = Path.of(outfile).getParent()
            assert ctx.readableDirs.contains(outParent)

            return dispatchResult
        } as AgentRunner

        when:
        final runner = new ScriptRunner([process: [executor: 'local'], workDir: work.toString()])
        runner.setScript(new ScriptFile(main))
        runner.execute()

        then:
        noExceptionThrown()
        and:
        captured != null
        dispatchResult != null
    }

    private static DispatchContext getDispatchContext() {
        // read the per-thread dispatch context the bridge sets before runner.run; the testRunner
        // closure executes on the SAME thread, so the context is visible there
        final f = ModuleToolBridge.getDeclaredField('CONTEXT')
        f.setAccessible(true)
        final ThreadLocal<DispatchContext> tl = (ThreadLocal<DispatchContext>) f.get(null)
        return tl.get()
    }

    def 'registry fetch failure falls back to meta.yml spec, still a per-module tool with enforced schema'() {
        given: 'a temp dir with a registry-installed module layout (has .module-info marker) + meta.yml'
        final dir = Files.createTempDirectory('test')
        final work = dir.resolve('work'); Files.createDirectories(work)

        // pre-install a module in the registry layout so recoverModuleRef succeeds
        final moduleDir = dir.resolve('modules').resolve('offline-test').resolve('mymod')
        Files.createDirectories(moduleDir)
        // write the .module-info marker so recoverModuleRef recognises it as a registry install
        moduleDir.resolve('.module-info').text = 'checksum=dummy'

        // classic untyped module with meta.yml (the fallback schema source)
        moduleDir.resolve('main.nf').text = '''\
            process MYMOD {
                input:
                tuple val(meta), path(fastq)

                output:
                tuple val(meta), path('out.txt')

                script:
                """
                touch out.txt
                """
            }
            '''.stripIndent()

        moduleDir.resolve('meta.yml').text = '''\
            name: MYMOD
            description: Offline test module
            input:
              - - name: meta
                  type: map
                  description: sample meta
                - name: fastq
                  type: file
                  description: input reads
            output:
              - - name: meta
                  type: map
                  description: sample meta
                - name: outfile
                  type: file
                  description: output file
            '''.stripIndent()

        final main = dir.resolve('main.nf')
        main.text = """\
            include { MYMOD } from './modules/offline-test/mymod/main.nf'

            agent assistant {
                model 'm'
                instruction 'i'
                tools 'module_run'

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
                assistant(channel.of('hi')).view { it }
            }
            """.stripIndent()

        and:
        AgentRunnerRequest captured = null
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req ->
            captured = req
            // per-module tool NAMED after the module even though registry fetch fails
            assert req.toolSpecs.size() == 1
            assert req.toolSpecs[0].name == 'MYMOD'
            assert req.toolSpecs*.name.every { it != 'module_run' }
            // the meta.yml fallback still yields an ENFORCED flattened schema
            final schema = req.toolSpecs[0].inputSchema
            assert (schema.properties as Map).containsKey('fastq')
            assert (schema.properties as Map).containsKey('meta')
            assert (schema.required as List).contains('fastq')
            assert schema.additionalProperties == false
            return 'done'
        } as AgentRunner

        when: 'the registry URL is unreachable so fetchModuleMetadata throws; fallback to meta.yml must happen'
        // configure an unreachable registry URL to force the fetch to throw an exception
        final runner = new ScriptRunner([
            process: [executor: 'local'],
            workDir: work.toString(),
            registry: [url: 'http://localhost:0']   // guaranteed to fail immediately
        ])
        runner.setScript(new ScriptFile(main))
        runner.execute()

        then: 'bridge is built using meta.yml spec without throwing'
        noExceptionThrown()
        and:
        captured != null
    }

    def 'module_run and filesystem capabilities together produce a per-module tool plus filesystem'() {
        given:
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req ->
            // the per-module tool (greet) + filesystem
            assert req.toolSpecs.size() == 2
            assert req.toolSpecs*.name.contains('greet')
            assert req.toolSpecs*.name.contains('filesystem')
            assert req.toolSpecs*.name.every { it != 'module_run' }
            // invoke the per-module tool to verify dispatch still works
            final result = req.dispatch.call('greet', '{"name":"World"}')
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
