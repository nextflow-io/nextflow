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

import nextflow.script.ScriptMeta
import spock.lang.Timeout
import test.Dsl2Spec

import static test.ScriptHelper.runScript

/**
 * Characterization tests for an {@code agent} authored as a MODULE and consumed with the
 * ordinary include statement, exactly like a process (design "agent module include" §3, §4.1).
 *
 * <p>The runtime include path ({@code IncludeDef.load0} -> {@code ScriptMeta.addModule} ->
 * {@code BaseScript.getComponent} -> {@code BindableDef.invoke_a}) is generic over
 * {@code ComponentDef} and never casts to {@code ProcessDef}, so an {@code AgentDef} rides it
 * unchanged. That behaviour was entirely unprotected by tests; this class locks it down.
 *
 * <p>The LLM is stubbed through {@link AgentRunnerProvider#testRunner} with a capturing runner,
 * so every assertion is deterministic and NO model call is made.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Timeout(120)
class AgentModuleIncludeTest extends Dsl2Spec {

    /** An agent module: no `nextflow.enable.types` needed for a typed `agent` block. */
    private static String agentModule(String name) {
        return """\
            agent ${name} {
                model 'openai/gpt-4o'
                instruction 'You write QA reports.'

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

    private List<AgentRunnerRequest> captured

    def setup() {
        captured = Collections.synchronizedList(new ArrayList<AgentRunnerRequest>())
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req ->
            captured.add(req)
            return "REPORT-${req.agentName}".toString()
        } as AgentRunner
    }

    def cleanup() {
        AgentRunnerProvider.testRunner = null
    }

    private static Path write(Path path, String text) {
        Files.createDirectories(path.parent)
        path.text = text
        return path
    }

    def 'should include an agent from a module file'() {
        given:
        final root = Files.createTempDirectory('test')
        write(root.resolve('mods/reporter.nf'), agentModule('reporter'))
        final main = write(root.resolve('main.nf'), '''\
            include { reporter } from './mods/reporter.nf'

            workflow {
                reporter(channel.of('s1')).view()
            }
            '''.stripIndent())

        when:
        final result = runScript(main)

        then:
        captured.size() == 1
        captured[0].agentName == 'reporter'
        result.val == 'REPORT-reporter'
    }

    // -- design §7.1/OQ-1: unlike a process, an agent has NO duplicate-invocation guard, and that
    //    was a DELIBERATE rejection (§15) -- AgentDef builds a fresh ProcessConfigV2 per call, so a
    //    second invocation binds its own channels and nothing is corrupted. Without this test the
    //    rejection lives only in a spec file, and a later refactor lifting the
    //    `this instanceof ProcessDef` gate in BindableDef.invoke_a up to BindableDef would silently
    //    turn a supported pattern into a DuplicateProcessInvocation error.
    def 'should allow two invocations of one included agent in the same scope'() {
        given:
        final root = Files.createTempDirectory('test')
        write(root.resolve('mods/reporter.nf'), agentModule('reporter'))
        final main = write(root.resolve('main.nf'), '''\
            include { reporter } from './mods/reporter.nf'

            workflow {
                reporter(channel.of('s1')).view()
                reporter(channel.of('s2')).view()
            }
            '''.stripIndent())

        when:
        runScript(main)

        then: 'both calls run -- no DuplicateProcessInvocation'
        noExceptionThrown()
        captured.size() == 2
        captured.collect { it.agentName } == ['reporter', 'reporter']
    }

    def 'should include an agent from a module directory via main.nf'() {
        given:
        final root = Files.createTempDirectory('test')
        write(root.resolve('mods/reporter/main.nf'), agentModule('reporter'))
        final main = write(root.resolve('main.nf'), '''\
            include { reporter } from './mods/reporter'

            workflow {
                reporter(channel.of('s1')).view()
            }
            '''.stripIndent())

        when:
        final result = runScript(main)

        then:
        captured.size() == 1
        captured[0].agentName == 'reporter'
        result.val == 'REPORT-reporter'
    }

    def 'should include an agent through a parent-relative path'() {
        given:
        final root = Files.createTempDirectory('test')
        write(root.resolve('mods/reporter/main.nf'), agentModule('reporter'))
        final main = write(root.resolve('pipeline/main.nf'), '''\
            include { reporter } from '../mods/reporter'

            workflow {
                reporter(channel.of('s1')).view()
            }
            '''.stripIndent())

        when:
        final result = runScript(main)

        then:
        captured[0].agentName == 'reporter'
        result.val == 'REPORT-reporter'
    }

    def 'should include an agent through an absolute path'() {
        given:
        final root = Files.createTempDirectory('test')
        final module = write(root.resolve('mods/reporter/main.nf'), agentModule('reporter'))
        final main = write(root.resolve('main.nf'), """\
            include { reporter } from '${module.toAbsolutePath()}'

            workflow {
                reporter(channel.of('s1')).view()
            }
            """.stripIndent())

        when:
        final result = runScript(main)

        then:
        captured[0].agentName == 'reporter'
        result.val == 'REPORT-reporter'
    }

    // -- ALIAS: cloneWithName renames the component, so the task/processor -- and therefore the
    //    work dir, the progress table and the trace -- carry the ALIAS, not the declared name.
    def 'should name the task after the alias when the include is aliased'() {
        given:
        final root = Files.createTempDirectory('test')
        write(root.resolve('mods/reporter/main.nf'), agentModule('reporter'))
        final main = write(root.resolve('main.nf'), '''\
            include { reporter as qc } from './mods/reporter'

            workflow {
                qc(channel.of('s1')).view()
            }
            '''.stripIndent())

        when:
        final result = runScript(main)

        then:
        captured[0].agentName == 'qc'
        result.val == 'REPORT-qc'
        and: 'both the declared name and the alias reach the config-selector registry, so a'
        // `withName:` block targeting either is not reported as unmatched by Session#checkConfig
        ScriptMeta.allAgentNames().contains('reporter')
        ScriptMeta.allAgentNames().contains('qc')
    }

    def 'should include two agents selectively from one module, one of them aliased'() {
        given:
        final root = Files.createTempDirectory('test')
        write(root.resolve('mods/agents.nf'), agentModule('alpha') + '\n' + agentModule('beta'))
        final main = write(root.resolve('main.nf'), '''\
            include { alpha ; beta as gamma } from './mods/agents.nf'

            workflow {
                alpha(channel.of('s1')).view()
                gamma(channel.of('s2')).view()
            }
            '''.stripIndent())

        when:
        runScript(main)

        then:
        (captured*.agentName as Set) == ['alpha', 'gamma'] as Set
    }

    // -- a module including ANOTHER module's agent resolves the path against its OWN location,
    //    and an invocation inside a named workflow qualifies the name with the workflow scope.
    def 'should resolve a nested module-to-module include and qualify the name inside a named workflow'() {
        given:
        final root = Files.createTempDirectory('test')
        write(root.resolve('mods/reporter/main.nf'), agentModule('qa'))
        write(root.resolve('mods/wrapper/main.nf'), '''\
            include { qa } from '../reporter'

            workflow wrapped {
                take:
                samples

                main:
                qa(samples).view()
            }
            '''.stripIndent())
        final main = write(root.resolve('main.nf'), '''\
            include { wrapped } from './mods/wrapper'

            workflow {
                wrapped(channel.of('s1'))
            }
            '''.stripIndent())

        when:
        runScript(main)

        then: 'the agent task is scoped by the enclosing named workflow'
        captured.size() == 1
        captured[0].agentName == 'wrapped:qa'
    }

    def 'should report a compile error for an unknown included agent name'() {
        given:
        final root = Files.createTempDirectory('test')
        final module = write(root.resolve('mods/reporter/main.nf'), agentModule('reporter'))
        final main = write(root.resolve('main.nf'), '''\
            include { nope } from './mods/reporter'

            workflow {
                nope(channel.of('s1')).view()
            }
            '''.stripIndent())

        when:
        runScript(main)

        then:
        final e = thrown(Exception)
        // note: the reported path is the REAL path, so match on the module-relative tail
        rootMessages(e).contains("Included name 'nope' is not defined in module '")
        rootMessages(e).contains(module.subpath(module.nameCount - 3, module.nameCount).toString())
        and: 'the agent is never run'
        captured.isEmpty()
    }

    def 'should report a compile error when the same agent name is included twice'() {
        given:
        final root = Files.createTempDirectory('test')
        write(root.resolve('mods/a/main.nf'), agentModule('reporter'))
        write(root.resolve('mods/b/main.nf'), agentModule('reporter'))
        final main = write(root.resolve('main.nf'), '''\
            include { reporter } from './mods/a'
            include { reporter } from './mods/b'

            workflow {
                reporter(channel.of('s1')).view()
            }
            '''.stripIndent())

        when:
        runScript(main)

        then:
        final e = thrown(Exception)
        rootMessages(e).contains('`reporter` is already included')
        and:
        captured.isEmpty()
    }

    def 'should report an error when the module directory has no main.nf'() {
        given:
        final root = Files.createTempDirectory('test')
        Files.createDirectories(root.resolve('mods/reporter'))
        final main = write(root.resolve('main.nf'), '''\
            include { reporter } from './mods/reporter'

            workflow {
                reporter(channel.of('s1')).view()
            }
            '''.stripIndent())

        when:
        runScript(main)

        then:
        thrown(Exception)
    }

    /** Flatten an exception chain (plus any suppressed causes) into one searchable string. */
    private static String rootMessages(Throwable e) {
        final sb = new StringBuilder()
        Throwable t = e
        while( t != null ) {
            sb.append(t.message ?: '').append('\n')
            for( final s : t.getSuppressed() )
                sb.append(s.message ?: '').append('\n')
            t = t.cause === t ? null : t.cause
        }
        return sb.toString()
    }
}
