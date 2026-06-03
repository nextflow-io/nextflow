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

/**
 * End-to-end test of a spec-driven module tool (Phase 3.2) using a CLASSIC DSL2
 * tuple module + a sibling {@code meta.yml}, faithful to the nf-core shape.
 *
 * The agent declares {@code tools './mod.nf'}; the bridge finds the sibling
 * {@code meta.yml}, derives a FLATTENED tool input schema ({@code meta}→object,
 * {@code reads}→string), and at dispatch time reassembles the flattened JSON args
 * into the {@code [meta, path]} tuple the process expects. The module runs through
 * the REAL local executor (it {@code cat}s the staged input into {@code out.txt}),
 * proving real staging + execution; the tuple/file output is serialized back to JSON
 * with the output file as an ABSOLUTE PATH STRING (the opaque-path contract).
 *
 * The {@code @Timeout} fails if the tool input queues are not poisoned on completion.
 */
@Timeout(90)
class AgentModuleSpecToolTest extends Dsl2Spec {

    def cleanup() {
        AgentRunnerProvider.testRunner = null
    }

    def 'should run a classic DSL2 tuple module (meta.yml-driven) as an agent tool and terminate'() {
        given:
        final dir = Files.createTempDirectory('test')
        final work = dir.resolve('work'); Files.createDirectories(work)
        final reads = dir.resolve('reads.txt'); reads.text = 'hello'
        final readsAbs = reads.toAbsolutePath().toString()

        // a CLASSIC DSL2 module mimicking the nf-core tuple shape (NOT typed records);
        // the script body uses `cat`, which runs via the LocalExecutor shell (no container).
        // The output is a `.dat` (NOT a text-like extension) so it stays an opaque path handle.
        dir.resolve('mod.nf').text = '''
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

        // sibling meta.yml describing the tuple I/O
        dir.resolve('meta.yml').text = '''\
            name: echo_tool
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

        final main = dir.resolve('main.nf')
        main.text = '''
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
                a(channel.of('go')).view { it }
            }
            '''.stripIndent()

        and:
        AgentRunnerRequest captured = null
        String dispatchResult = null
        Map<String,Object> outAssert = [:]
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req ->
            captured = req
            // the bridge exposes the spec-driven `echo_tool` with a FLATTENED schema
            assert req.toolSpecs.size() == 1
            assert req.toolSpecs[0].name == 'echo_tool'
            assert req.toolSpecs[0].inputSchema.properties.meta.type == 'object'
            assert req.toolSpecs[0].inputSchema.properties.reads.type == 'string'
            // invoke the tool: this drives the REAL echo_tool process through the executor,
            // staging `reads.txt` and producing `out.txt`
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
            // the agent's final answer
            return dispatchResult
        } as AgentRunner

        when:
        final runner = new ScriptRunner([process: [executor: 'local'], workDir: work.toString()])
        runner.setScript(new ScriptFile(main))
        runner.execute()

        then:
        // the dispatch went through the bridge and returned the real process output
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

    def 'should inline a small structured json tool output to the LLM instead of a path handle'() {
        given:
        final dir = Files.createTempDirectory('test')
        final work = dir.resolve('work'); Files.createDirectories(work)
        final reads = dir.resolve('reads.txt'); reads.text = 'hello'
        final readsAbs = reads.toAbsolutePath().toString()

        // a CLASSIC DSL2 module producing a SMALL .json stats output (the assemblyscan
        // shape): the script body writes a tiny JSON that the LLM must reason over
        dir.resolve('mod.nf').text = '''
            process echo_tool {
                input:
                tuple val(meta), path(reads)

                output:
                tuple val(meta), path("stats.json"), emit: report

                script:
                """
                echo '{"n50":54321,"contigs":12}' > stats.json
                """
            }
            '''.stripIndent()

        // sibling meta.yml describing the tuple I/O; the file output is a .json
        dir.resolve('meta.yml').text = '''\
            name: echo_tool
            description: Compute assembly statistics
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
                  description: the stats json
                  pattern: "*.json"
            '''.stripIndent()

        final main = dir.resolve('main.nf')
        main.text = '''
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
                a(channel.of('go')).view { it }
            }
            '''.stripIndent()

        and:
        String dispatchResult = null
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req ->
            dispatchResult = req.dispatch.call('echo_tool', JsonOutput.toJson([meta: [id: 's1'], reads: readsAbs]))
            final parsed = new JsonSlurper().parseText(dispatchResult) as Map
            // the result is keyed by the emit name `report`; its `outfile` carries the
            // FILE CONTENTS as a String (inlined) -- NOT an absolute path handle
            assert parsed.containsKey('report')
            final report = parsed.report as Map
            assert report.meta == [id: 's1']
            final outfile = report.outfile as String
            assert !outfile.startsWith('/')
            assert outfile.contains('"n50"')
            assert outfile.contains('54321')
            return dispatchResult
        } as AgentRunner

        when:
        final runner = new ScriptRunner([process: [executor: 'local'], workDir: work.toString()])
        runner.setScript(new ScriptFile(main))
        runner.execute()

        then:
        dispatchResult != null
        final parsed = new JsonSlurper().parseText(dispatchResult) as Map
        parsed.containsKey('report')
        (parsed.report as Map).meta == [id: 's1']
        and:
        // the file contents were inlined: the value is the JSON content, not a /abs/path
        final outfile = (parsed.report as Map).outfile as String
        !outfile.startsWith('/')
        outfile.contains('"n50"')
        outfile.contains('54321')
    }

    def 'should skip an nf-core versions (eval/topic) output and not block the dispatch'() {
        given:
        final dir = Files.createTempDirectory('test')
        final work = dir.resolve('work'); Files.createDirectories(work)
        final reads = dir.resolve('reads.txt'); reads.text = 'hello'
        final readsAbs = reads.toAbsolutePath().toString()

        // a module with a DATA output (`report`) AND an nf-core `versions` output:
        // a tuple carrying an `eval` component routed to a `topic` -- exactly the
        // skesa shape that previously blocked the dispatcher forever on `.val`.
        // The data output is a `.dat` (NOT text-like) so it stays an opaque path handle.
        dir.resolve('mod.nf').text = '''
            process echo_tool {
                input:
                tuple val(meta), path(reads)

                output:
                tuple val(meta), path("out.dat"), emit: report
                tuple val("${task.process}"), val('echo'), eval('echo 1.0'), topic: versions, emit: versions_echo

                script:
                """
                cat ${reads} > out.dat
                """
            }
            '''.stripIndent()

        // sibling meta.yml: a `report` data output + a versions output whose components
        // include an `eval` (plus the matching `versions` topic)
        dir.resolve('meta.yml').text = '''\
            name: echo_tool
            description: Copy reads to an output file
            input:
              - - name: meta
                  type: map
                - name: reads
                  type: file
            output:
              - - name: meta
                  type: map
                - name: outfile
                  type: file
              - - name: proc
                  type: string
                - name: tool
                  type: string
                - name: version
                  type: eval
            topics:
              - - name: proc
                  type: string
                - name: tool
                  type: string
                - name: version
                  type: eval
            '''.stripIndent()

        dir.resolve('main.nf').text = '''
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
                a(channel.of('go')).view { it }
            }
            '''.stripIndent()

        and:
        String dispatchResult = null
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req ->
            dispatchResult = req.dispatch.call('echo_tool', JsonOutput.toJson([meta: [id: 's1'], reads: readsAbs]))
            return dispatchResult
        } as AgentRunner

        when:
        final runner = new ScriptRunner([process: [executor: 'local'], workDir: work.toString()])
        runner.setScript(new ScriptFile(dir.resolve('main.nf')))
        runner.execute()

        then:
        // the dispatch returned (it did NOT block forever on the eval/versions output's `.val`)
        // and the whole run terminated within @Timeout
        dispatchResult != null
        final parsed = new JsonSlurper().parseText(dispatchResult) as Map
        // the data output IS collected ...
        parsed.containsKey('report')
        ((parsed.report as Map).outfile as String).startsWith('/')
        // ... and the eval/versions output is NOT (skipped, not blocked on)
        !parsed.containsKey('versions_echo')
        parsed.size() == 1
    }
}
