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
import spock.lang.Timeout
import test.Dsl2Spec

import static test.ScriptHelper.runScript

/**
 * The tool scope of a MODULE agent is MODULE-LEXICAL (design §5): an agent sees only the
 * processes defined or included by the script that DECLARES it. The including script's processes
 * are not visible, exactly as a process body can only call what its own script defines or
 * includes. That keeps the LLM's tool surface -- and therefore the model's behaviour -- independent
 * of who imported the agent, which is what makes {@code main.nf} + {@code skills/} + {@code tools/}
 * a shippable unit.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Timeout(120)
class AgentModuleRunScopeTest extends Dsl2Spec {

    private AgentRunnerRequest captured

    def cleanup() {
        AgentRunnerProvider.testRunner = null
    }

    private static Path write(Path path, String text) {
        Files.createDirectories(path.parent)
        path.text = text
        return path
    }

    private static String agentModule(String toolsDirective, String extra = '') {
        return """\
            ${extra}
            agent reporter {
                model 'openai/gpt-4o'
                instruction 'You write QA reports.'
                ${toolsDirective}

                input:
                sample: String

                output:
                report: String

                prompt:
                \"\"\"
                Report on \${sample}.
                \"\"\"
            }
            """.stripIndent()
    }

    /** A typed `exec` process module; `nextflow.enable.types` is per-FILE, hence the flag. */
    private static String toolModule(String name, String prefix) {
        return """\
            nextflow.enable.types = true

            process ${name} {
                input:
                text: String

                output:
                out: String

                exec:
                out = '${prefix}' + text
            }
            """.stripIndent()
    }

    def 'nf:module_run advertises the process the MODULE includes, not one defined only in the caller'() {
        given:
        final root = Files.createTempDirectory('test')
        write(root.resolve('mods/reporter/tools/shouter.nf'), toolModule('shouter', 'MOD:'))
        write(root.resolve('mods/reporter/main.nf'), agentModule(
            "tools 'nf:module_run'",
            "include { shouter } from './tools/shouter.nf'\n"))
        final main = write(root.resolve('main.nf'), '''\
            nextflow.enable.types = true

            include { reporter } from './mods/reporter'

            process caller_only {
                input:
                text: String

                output:
                out: String

                exec:
                out = 'CALLER:' + text
            }

            workflow {
                reporter(channel.of('s1')).view()
            }
            '''.stripIndent())

        and:
        String dispatched = null
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req ->
            captured = req
            dispatched = req.dispatch.call('shouter', '{"text":"hi"}')
            return 'ok'
        } as AgentRunner

        when:
        runScript(main)

        then: 'the module-included process IS advertised'
        captured.toolSpecs*.name == ['shouter']
        and: 'the process defined only in the INCLUDING script is NOT'
        !captured.toolSpecs*.name.contains('caller_only')
        and: 'a full round trip through the module`s own tool works'
        new JsonSlurper().parseText(dispatched) == [out: 'MOD:hi']
    }

    def 'a process included only by the CALLER is not a member of the module agent`s nf:module_run'() {
        given: 'the module includes `whisper`, the entry script includes `shouter`'
        final root = Files.createTempDirectory('test')
        write(root.resolve('mods/reporter/tools/whisper.nf'), toolModule('whisper', 'MOD:'))
        write(root.resolve('mods/reporter/main.nf'), agentModule(
            "tools 'nf:module_run:shouter'",
            "include { whisper } from './tools/whisper.nf'\n"))
        write(root.resolve('tools/shouter.nf'), toolModule('shouter', 'CALLER:'))
        final main = write(root.resolve('main.nf'), '''\
            include { reporter } from './mods/reporter'
            include { shouter } from './tools/shouter.nf'

            workflow {
                reporter(channel.of('s1')).view()
            }
            '''.stripIndent())

        and:
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req -> 'ok' } as AgentRunner

        when: 'the named tool is included ONLY by the entry script'
        runScript(main)

        then: 'G8(c) - an explicit leaf that does not exist is a hard error'
        final e = thrown(Exception)
        final msg = messages(e)
        msg.contains('Agent `reporter`: Tool `nf:module_run:shouter` does not exist')
        and: 'the available list is the MODULE`s scope alone, which is what makes the rule visible'
        msg.contains('available: `nf:module_run:whisper`')
    }

    def 'a module agent relative include resolves against the module directory'() {
        given:
        final root = Files.createTempDirectory('test')
        write(root.resolve('mods/reporter/tools/wordcount.nf'), toolModule('wordcount', 'MOD:'))
        and: 'a DECOY of the same relative path beside the entry script'
        write(root.resolve('tools/wordcount.nf'), toolModule('wordcount', 'DECOY:'))
        write(root.resolve('mods/reporter/main.nf'), agentModule(
            "tools 'nf:module_run:wordcount'",
            "include { wordcount } from './tools/wordcount.nf'\n"))
        final main = write(root.resolve('main.nf'), '''\
            include { reporter } from './mods/reporter'

            workflow {
                reporter(channel.of('s1')).view()
            }
            '''.stripIndent())

        and:
        String dispatched = null
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req ->
            captured = req
            dispatched = req.dispatch.call('wordcount', '{"text":"hi"}')
            return 'ok'
        } as AgentRunner

        when:
        runScript(main)

        then:
        captured.toolSpecs*.name == ['wordcount']
        and: 'the MODULE`s copy ran, not the decoy beside the entry script'
        new JsonSlurper().parseText(dispatched) == [out: 'MOD:hi']
    }

    private static String messages(Throwable e) {
        final sb = new StringBuilder()
        Throwable t = e
        while( t != null ) {
            sb.append(t.message ?: '').append('\n')
            t = t.cause === t ? null : t.cause
        }
        return sb.toString()
    }
}
