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

import groovy.json.JsonOutput
import groovy.json.JsonSlurper
import nextflow.script.ScriptFile
import nextflow.script.ScriptRunner
import spock.lang.Timeout
import test.Dsl2Spec

/**
 * Regression test for the topic-routed {@code versions} output bug: an nf-core module whose
 * registry {@code meta.yml} types the eval/version component as {@code string} (NOT {@code eval},
 * e.g. nf-core/assemblyscan) used to BLOCK the agent dispatch forever — the bridge read the
 * topic-source channel's {@code .val}, which never binds a per-invocation value. The fix detects
 * topic-routed outputs from the ProcessDef ({@code OutParam.getChannelTopicName()}), which is
 * authoritative; the meta.yml {@code type} is unreliable.
 *
 * <p>This lives in its OWN spec class on purpose: running two {@code topic}-creating agent sessions
 * back-to-back in the same JVM trips a pre-existing GPars cross-session operator-join hang at
 * session shutdown (run agent test classes individually). Keeping at most one topic-agent-session
 * per class avoids that — {@link AgentModuleSpecToolTest} already carries the {@code eval}-typed
 * topic case.
 */
@Timeout(90)
class AgentTopicOutputToolTest extends Dsl2Spec {

    def cleanup() {
        AgentRunnerProvider.testRunner = null
    }

    def 'should skip a topic-routed versions output typed as string and not block the dispatch'() {
        given:
        final dir = Files.createTempDirectory('test')
        final work = dir.resolve('work'); Files.createDirectories(work)
        final reads = dir.resolve('reads.txt'); reads.text = 'hello'
        final readsAbs = reads.toAbsolutePath().toString()

        // The ASSEMBLYSCAN shape: a DATA output (`report`, a small .json that is INLINED) AND
        // an nf-core `versions` output routed to a `topic`. Unlike skesa, assemblyscan's
        // registry-generated meta.yml types the version-command component as `type: string`
        // (NOT `eval`), so the old `isEvalOutput` MISSES it -- the dispatcher used to block
        // forever reading the topic-source channel's `.val`. The fix skips it via the
        // ProcessDef's topic-routed output param, which is authoritative.
        dir.resolve('mod.nf').text = '''
            process echo_tool {
                input:
                tuple val(meta), path(reads)

                output:
                tuple val(meta), path("out.json"), emit: report
                tuple val("${task.process}"), val('mytool'), eval('echo 1.0'), topic: versions, emit: versions_mytool

                script:
                """
                echo '{"n50":42}' > out.json
                """
            }
            '''.stripIndent()

        // sibling meta.yml: a `report` data output + a versions output whose components are
        // ALL typed `string` (NO `eval` anywhere) -- mirroring assemblyscan's real meta.yml
        dir.resolve('meta.yml').text = '''\
            name: echo_tool
            description: Compute assembly statistics
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
                  pattern: "*.json"
              - - name: proc
                  type: string
                - name: tool
                  type: string
                - name: version
                  type: string
            topics:
              - - name: proc
                  type: string
                - name: tool
                  type: string
                - name: version
                  type: string
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
        // the dispatch returned (it did NOT block forever on the string-typed topic output's
        // `.val`) and the whole run terminated within @Timeout
        dispatchResult != null
        final parsed = new JsonSlurper().parseText(dispatchResult) as Map
        // the data output IS collected ...
        parsed.containsKey('report')
        and:
        // ... and it was INLINED (small .json): the `outfile` carries the JSON CONTENTS
        final report = parsed.report as Map
        report.meta == [id: 's1']
        final outfile = report.outfile as String
        !outfile.startsWith('/')
        outfile.contains('n50')
        outfile.contains('42')
        and:
        // ... and the topic-routed versions output is NOT present (skipped, not blocked on)
        !parsed.containsKey('versions_mytool')
        parsed.size() == 1
    }
}
