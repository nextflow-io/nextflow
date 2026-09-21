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

import nextflow.exception.ScriptRuntimeException
import nextflow.script.AgentDef
import spock.lang.Timeout
import spock.lang.Unroll
import test.Dsl2Spec

import static test.ScriptHelper.runScript

/**
 * The two checks that make the tool namespace safe rather than merely expressive (§4), plus the
 * runner-native contribution to the resume key (§7).
 *
 * <p>The declaration grammar guarantees that a ref an author WRITES can always become a wire name
 * — {@link ToolRef} admits only {@code [A-Za-z0-9_-]} in a segment. What it cannot guarantee is
 * what a glob or a bare non-leaf DRAGS IN: {@code nf:module_run} enumerates every in-scope
 * process, and a Nextflow process name is a {@code JavaLetter JavaLetterOrDigit*}, which admits
 * {@code $}, non-ASCII letters and any length. Those are the names checked here.
 *
 * <p>The collision half covers the namespace the MODEL sees, which is wider than the declared one:
 * the runner injects {@code activate_skill}/{@code read_skill_resource} for a skills agent and
 * {@code final_answer} for a structured-output agent on the canonical runner, and the model cannot
 * tell any of them apart from a tool.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Timeout(60)
class ToolNameValidationTest extends Dsl2Spec {

    def cleanup() {
        AgentRunnerProvider.testRunner = null
    }

    /** The selection exactly as {@code AgentDef} builds it: the grammar over the in-scope processes. */
    private static ToolRefResolver.Selection select(List<String> processNames, List<String> declared) {
        return ToolRefResolver.standard('Agent `assistant`', processNames).resolve(declared)
    }

    /**
     * The §4 pass over one selection. The named options are the state of the AGENT that decides
     * which tools the runner injects beside the declared ones.
     */
    private static void check(Map opts = [:], ToolRefResolver.Selection selection) {
        AgentDef.checkWireNames('assistant', selection,
                opts.skills as boolean, opts.structuredOutput as boolean, opts.containerized as boolean)
    }

    // -----------------------------------------------------------------------
    // §4 — validate, never sanitize
    // -----------------------------------------------------------------------

    def 'should reject a selected process whose name the wire namespace cannot carry'() {
        given: 'a process name that is legal Nextflow and illegal as an OpenAI function name'
        def selection = select(['my$proc'], ['nf:module_run'])

        when:
        check(selection)

        then:
        def e = thrown(ScriptRuntimeException)
        and: 'the process is named, and so is the character that made it illegal'
        e.message.contains('my$proc')
        e.message.contains('the illegal character(s) `$`')
        and: 'the rule the name has to satisfy is quoted'
        e.message.contains('[a-zA-Z0-9_-]')
        and: 'the ref it was dragged in by is named, so the author knows which entry to narrow'
        e.message.contains('nf:module_run:my$proc')
        and: 'a rename is refused, not performed'
        e.message.contains('never rewritten automatically')
    }

    def 'should never silently rename an illegal process into an existing one'() {
        given: 'the collision a sanitizing implementation would create: my$proc -> my_proc'
        def selection = select(['my_proc', 'my$proc'], ['nf:module_run'])

        when:
        check(selection)

        then: 'the illegal name is an error; it is never merged into its legal neighbour'
        def e = thrown(ScriptRuntimeException)
        e.message.contains('my$proc')
        and: 'and the selection is NOT quietly reduced to the legal one'
        selection.brokeredNames == ['my_proc', 'my$proc']
    }

    def 'should reject a selected process whose name exceeds the wire length limit'() {
        when:
        check(select(['P' * 65], ['nf:module_run:*']))

        then: 'the limit is quoted together with the length that broke it'
        def e = thrown(ScriptRuntimeException)
        e.message.contains('65 characters (the limit is 64)')

        when: 'one character shorter'
        check(select(['P' * 64], ['nf:module_run:*']))

        then: 'the boundary is inclusive'
        noExceptionThrown()
    }

    @Unroll
    def 'should reject the illegal wire name #NAME'() {
        when:
        check(select([NAME], ['nf:module_run']))

        then:
        def e = thrown(ScriptRuntimeException)
        e.message.contains(NAME)
        e.message.contains(OFFENDER)

        where:
        NAME        | OFFENDER
        'my$proc'   | '`$`'
        'Σ_SORT'    | '`Σ`'
        'a.b'       | '`.`'
        'a b'       | '` `'
        'naïve$sum' | '`ï`'      // every offending character is listed, not just the first
    }

    @Unroll
    def 'should accept the legal wire name #NAME'() {
        when:
        check(select([NAME], ['nf:module_run']))

        then:
        noExceptionThrown()

        where:
        NAME << ['GREET', 'greet', 'SAMTOOLS_SORT', 'a-b', 'x9', '_leading', 'P' * 64]
    }

    def 'should validate the runner-native names too, not only the processes'() {
        when: 'the fixed fs:/shell: leaves go through the same check'
        check(select([], ['fs:*']))

        then: 'they are legal by construction — this is the regression guard on that claim'
        noExceptionThrown()
    }

    // -----------------------------------------------------------------------
    // §4 — one wire name, one source
    // -----------------------------------------------------------------------

    def 'should reject a process colliding with a runner-native tool, naming both sources'() {
        when:
        check(select(['read', 'greet'], ['nf:module_run:read', 'fs:read']))

        then:
        def e = thrown(ScriptRuntimeException)
        e.message.contains('`read`')
        and: 'BOTH sources are named — the message says which two things claim the name'
        e.message.contains('`fs:read`')
        e.message.contains('`nf:module_run:read`')
        and: 'the reason the flat namespace forbids it'
        e.message.contains('single flat namespace')
    }

    def 'should report a collision identically whatever order the refs were declared in'() {
        when:
        check(select(['read'], ['nf:module_run:read', 'fs:read']))
        then:
        def first = thrown(ScriptRuntimeException)

        when: 'the very same directive, written the other way round'
        check(select(['read'], ['fs:read', 'nf:module_run:read']))
        then:
        def second = thrown(ScriptRuntimeException)

        and: 'byte-identical — the sources are sorted, never listed in declaration order'
        first.message == second.message
        and: 'sorted, so `fs:read` precedes `nf:module_run:read`'
        first.message.indexOf('`fs:read`') < first.message.indexOf('`nf:module_run:read`')
    }

    def 'should not report a collision when a glob and an explicit ref select the same tool'() {
        when: 'G9: an identical (wire name, source) pair collapses silently'
        check(select(['GREET'], ['nf:module_run', 'nf:module_run:GREET']))
        then:
        noExceptionThrown()

        when:
        check(select([], ['fs:*', 'fs:read']))
        then:
        noExceptionThrown()
    }

    def 'should reject a process colliding with the skills tools the runner injects'() {
        given:
        def selection = select(['activate_skill', 'greet'], ['nf:module_run'])

        when: 'the agent declares no skills, so the runner injects nothing'
        check(selection)
        then: 'the name is unremarkable — it is a tool like any other'
        noExceptionThrown()

        when: 'the same agent declares skills'
        check(selection, skills: true)
        then:
        def e = thrown(ScriptRuntimeException)
        e.message.contains('`activate_skill`')
        e.message.contains('`nf:module_run:activate_skill`')
        e.message.contains('the `skills` directive')
    }

    def 'should reject a process colliding with read_skill_resource'() {
        when:
        check(select(['read_skill_resource'], ['nf:module_run']), skills: true)

        then:
        def e = thrown(ScriptRuntimeException)
        e.message.contains('`read_skill_resource`')
        e.message.contains('the `skills` directive')
    }

    def 'should reject a process named final_answer on the canonical runner'() {
        given: 'the tool the canonical runner injects to terminate a structured-output agent'
        def selection = select(['final_answer'], ['nf:module_run'])

        when:
        check(selection, structuredOutput: true, containerized: true)

        then:
        def e = thrown(ScriptRuntimeException)
        e.message.contains('`final_answer`')
        e.message.contains('`nf:module_run:final_answer`')
        e.message.contains('the agent output declaration')

        when: 'the same agent with a free-text output declares no such tool'
        check(selection, containerized: true)
        then:
        noExceptionThrown()

        when: 'and the in-JVM runner decodes the structured answer without injecting a tool'
        check(selection, structuredOutput: true)
        then: 'so the name is free there — this is a per-runner namespace, like shell:bash'
        noExceptionThrown()
    }

    def 'should tolerate the injected names when nothing else claims them'() {
        when: 'the two skills tools share a source, so they never collide with each other'
        check(select(['greet'], ['nf:module_run']), skills: true, structuredOutput: true, containerized: true)

        then:
        noExceptionThrown()
    }

    def 'should check the injected names even for an agent with no tools at all'() {
        when:
        check(null, skills: true, structuredOutput: true, containerized: true)

        then:
        noExceptionThrown()
    }

    // -----------------------------------------------------------------------
    // §4 — wired into the real agent-build path, through the DSL
    // -----------------------------------------------------------------------

    def 'a process colliding with an fs: tool fails the agent at build time'() {
        given:
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req -> 'never reached' } as AgentRunner

        when:
        runScript('''
            nextflow.enable.types = true

            process read {
                input:
                name: String

                output:
                text: String

                exec:
                text = "read ${name}"
            }

            agent assistant {
                model 'm'
                instruction 'i'
                tools 'nf:module_run:read', 'fs:read'

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
        def e = thrown(Exception)
        def msg = messages(e)
        msg.contains('the tool name `read` is claimed by')
        msg.contains('`fs:read`')
        msg.contains('`nf:module_run:read`')
    }

    def 'a process whose name is illegal on the wire fails the agent at build time'() {
        given:
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req -> 'never reached' } as AgentRunner
        and: 'a process name that is legal Nextflow but 65 characters long'
        def long65 = 'P' * 65

        when: 'it is dragged in by a bare non-leaf ref rather than named explicitly'
        runScript('''
            nextflow.enable.types = true

            process LONG_NAME {
                input:
                name: String

                output:
                text: String

                exec:
                text = name
            }

            agent assistant {
                model 'm'
                instruction 'i'
                tools 'nf:module_run'

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
            '''.replace('LONG_NAME', long65))

        then:
        def e = thrown(Exception)
        def msg = messages(e)
        msg.contains('65 characters (the limit is 64)')
        msg.contains('never rewritten automatically')
    }

    // -----------------------------------------------------------------------
    // §7 — the runner-native contribution to the resume key
    // -----------------------------------------------------------------------

    def 'a runner-native tool must change the cache key, having no descriptor to hash'() {
        when: 'an agent whose only tools are runner-native still fingerprints'
        def withFs = AgentDef.toolsFingerprint(null, null, ['fs:read'], 'pi@1.0.0')

        then: 'the §7 hole: on the descriptor-only form this would have been null'
        withFs != null
        and: 'and it is NOT the key of the same agent with no tools'
        withFs != AgentDef.toolsFingerprint(null, null, null, 'pi@1.0.0')
    }

    def 'two agents differing only in their native tools get different keys'() {
        given:
        def runner = 'pi@1.0.0'

        expect: 'gaining a tool re-runs the agent rather than replaying a stale generation'
        AgentDef.toolsFingerprint(null, null, ['fs:read'], runner) !=
            AgentDef.toolsFingerprint(null, null, ['fs:read', 'fs:write'], runner)
        and: 'a different tool is a different key'
        AgentDef.toolsFingerprint(null, null, ['fs:read'], runner) !=
            AgentDef.toolsFingerprint(null, null, ['fs:write'], runner)
        and: 'declaring shell:bash is visible in the key'
        AgentDef.toolsFingerprint(null, null, ['fs:read'], runner) !=
            AgentDef.toolsFingerprint(null, null, ['fs:read', 'shell:bash'], runner)
        and: 'the same set is the same key, whatever order it was resolved in'
        AgentDef.toolsFingerprint(null, null, ['fs:read', 'fs:write'], runner) ==
            AgentDef.toolsFingerprint(null, null, ['fs:write', 'fs:read'], runner)
    }

    def 'the same agent gets a different key on a different runner'() {
        expect: 'the same fs:read is a different implementation on each runner'
        AgentDef.toolsFingerprint(null, null, ['fs:read'], 'pi@1.0.0') !=
            AgentDef.toolsFingerprint(null, null, ['fs:read'], 'langchain4j@1.0.0')

        and: 'and a runner UPGRADE invalidates it too — the version pins the image, which pins the SDK'
        AgentDef.toolsFingerprint(null, null, ['fs:read'], 'pi@1.0.0') !=
            AgentDef.toolsFingerprint(null, null, ['fs:read'], 'pi@1.1.0')

        and: 'an unchanged runner keeps the key stable'
        AgentDef.toolsFingerprint(null, null, ['fs:read'], 'pi@1.0.0') ==
            AgentDef.toolsFingerprint(null, null, ['fs:read'], 'pi@1.0.0')
    }

    def 'the runner identity is folded in only when a native tool needs it'() {
        given: 'a brokered-only agent — its tools are Nextflow tasks, identical on every runner'
        def brokered = [new ToolDescriptor('greet', 'say hi', [type: 'object'], null)]

        expect: 'the runner argument does not move a brokered-only key'
        AgentDef.toolsFingerprint(brokered, [greet: 'body'], null, 'pi@1.0.0') ==
            AgentDef.toolsFingerprint(brokered, [greet: 'body'], null, 'langchain4j@9.9.9')
        and: 'which is exactly the 2-arg form, byte for byte — no existing agent`s key moves'
        AgentDef.toolsFingerprint(brokered, [greet: 'body'], null, null) ==
            AgentDef.toolsFingerprint(brokered, [greet: 'body'])
        and: 'a tool-free agent still fingerprints to null, keeping the canonical source unchanged'
        AgentDef.toolsFingerprint(null, null, null, 'pi@1.0.0') == null
        AgentDef.toolsFingerprint([] as List<ToolDescriptor>, [:], [] as List<String>, 'pi@1.0.0') == null
    }

    def 'a mixed agent keys on both halves'() {
        given:
        def brokered = [new ToolDescriptor('greet', 'say hi', [type: 'object'], null)]
        def sources = [greet: 'body']
        def mixed = AgentDef.toolsFingerprint(brokered, sources, ['fs:read'], 'pi@1.0.0')

        expect: 'neither half alone is the same key'
        mixed != AgentDef.toolsFingerprint(brokered, sources, null, 'pi@1.0.0')
        mixed != AgentDef.toolsFingerprint(null, null, ['fs:read'], 'pi@1.0.0')
        and: 'and editing the brokered process still invalidates it'
        mixed != AgentDef.toolsFingerprint(brokered, [greet: 'edited'], ['fs:read'], 'pi@1.0.0')
    }

    def 'the runner identity degrades to the bare name when no plugin owns the runner'() {
        given: 'the AgentRunnerProvider test seam injects a runner with no plugin behind it'
        def runner = { AgentRunnerRequest req -> 'ok' } as AgentRunner

        when:
        def id = AgentDef.runnerIdentity(runner)

        then: 'a missing version is never made up, and never throws'
        id == runner.getName()
        and: 'it is deterministic within an installation, which is all the cache key needs'
        id == AgentDef.runnerIdentity(runner)
        and:
        AgentDef.runnerIdentity(null) == null
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
