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

import spock.lang.Timeout
import test.Dsl2Spec

import static test.ScriptHelper.runScript

/**
 * The KEY requirement of the agent-module work: a module agent's {@code skills} resolve from the
 * dir of the module that DECLARES the agent ({@code <moduleDir>/skills/<name>/SKILL.md}), never
 * from the including script's directory -- and that holds under an ALIAS too, because
 * {@code AgentDef.cloneWithName} is a shallow clone that preserves {@code owner}, and
 * {@code AgentDef.ownerBaseDir} prefers {@code ScriptMeta.getModuleDir()} over
 * {@code session.baseDir}.
 *
 * <p>Oracle note: skills do NOT travel in the prompt. Core renders the prompt from the
 * {@code PromptDef} alone; skills reach the runner only as {@code AgentRunnerRequest.skills}
 * (system-message injection happens in the nf-agent plugin), so the assertions here read
 * {@code captured.skills} from a capturing stub runner -- no model call.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Timeout(120)
class AgentModuleSkillsTest extends Dsl2Spec {

    private static final String MODULE_MARKER = 'MODULE-SKILL-BODY'
    private static final String STYLE_MARKER = 'MODULE-STYLE-BODY'
    private static final String DECOY_MARKER = 'DECOY-SKILL-BODY'

    private AgentRunnerRequest captured

    def setup() {
        captured = null
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req ->
            captured = req
            return 'ok'
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

    /** A `<dir>/SKILL.md` with YAML front matter, the shape SkillResolver.parseSkillDir expects. */
    private static void skill(Path skillDir, String name, String body) {
        write(skillDir.resolve('SKILL.md'), """\
            ---
            name: ${name}
            description: The ${name} skill
            ---

            ${body}
            """.stripIndent())
    }

    private static String agentModule(String name, String skillsDirective) {
        return """\
            agent ${name} {
                model 'openai/gpt-4o'
                instruction 'You write QA reports.'
                ${skillsDirective}

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

    private static String entryScript(String includeStmt, String call) {
        return """\
            ${includeStmt}

            workflow {
                ${call}(channel.of('s1')).view()
            }
            """.stripIndent()
    }

    def 'should resolve a module agent skill from the module directory'() {
        given:
        final root = Files.createTempDirectory('test')
        final mod = root.resolve('mods/reporter')
        write(mod.resolve('main.nf'), agentModule('reporter', "skills 'qa-report'"))
        skill(mod.resolve('skills/qa-report'), 'qa-report', MODULE_MARKER)
        final main = write(root.resolve('main.nf'),
            entryScript("include { reporter } from './mods/reporter'", 'reporter'))

        when:
        runScript(main)

        then:
        captured != null
        captured.skills*.name == ['qa-report']
        captured.skills[0].content.contains(MODULE_MARKER)
        and: 'the skill is NOT smuggled through the prompt'
        !captured.prompt.contains(MODULE_MARKER)
    }

    // -- the alias-safety case: cloneWithName preserves `owner`, so the skills root does not move
    def 'should resolve module skills from the module dir even when the agent is aliased'() {
        given:
        final root = Files.createTempDirectory('test')
        final mod = root.resolve('mods/reporter')
        write(mod.resolve('main.nf'), agentModule('reporter', "skills 'qa-report'"))
        skill(mod.resolve('skills/qa-report'), 'qa-report', MODULE_MARKER)
        final main = write(root.resolve('main.nf'),
            entryScript("include { reporter as qc } from './mods/reporter'", 'qc'))

        when:
        runScript(main)

        then:
        captured.agentName == 'qc'
        captured.skills*.name == ['qa-report']
        captured.skills[0].content.contains(MODULE_MARKER)
    }

    // -- "one or more skills", way 1: several entries under the same module skills root
    def 'should resolve two skills declared as two entries in the module directory'() {
        given:
        final root = Files.createTempDirectory('test')
        final mod = root.resolve('mods/reporter')
        write(mod.resolve('main.nf'), agentModule('reporter', "skills 'qa-report', 'style'"))
        skill(mod.resolve('skills/qa-report'), 'qa-report', MODULE_MARKER)
        skill(mod.resolve('skills/style'), 'style', STYLE_MARKER)
        final main = write(root.resolve('main.nf'),
            entryScript("include { reporter as qc } from './mods/reporter'", 'qc'))

        when:
        runScript(main)

        then:
        captured.skills*.name == ['qa-report', 'style']
        captured.skills.find { it.name == 'qa-report' }.content.contains(MODULE_MARKER)
        captured.skills.find { it.name == 'style' }.content.contains(STYLE_MARKER)
    }

    // -- "one or more skills", way 2: ONE entry whose dir has no SKILL.md expands to every
    //    immediate subdirectory that has one (SkillResolver.scanSkillRoot)
    def 'should expand a single skills entry whose directory holds several skills'() {
        given:
        final root = Files.createTempDirectory('test')
        final mod = root.resolve('mods/reporter')
        write(mod.resolve('main.nf'), agentModule('reporter', "skills 'bundle'"))
        skill(mod.resolve('skills/bundle/qa-report'), 'qa-report', MODULE_MARKER)
        skill(mod.resolve('skills/bundle/style'), 'style', STYLE_MARKER)
        final main = write(root.resolve('main.nf'),
            entryScript("include { reporter } from './mods/reporter'", 'reporter'))

        when:
        runScript(main)

        then:
        (captured.skills*.name as Set) == ['qa-report', 'style'] as Set
        captured.skills*.content.join('').contains(MODULE_MARKER)
        captured.skills*.content.join('').contains(STYLE_MARKER)
    }

    // -- THE negative case: a same-named `skills/` beside the ENTRY script must be IGNORED.
    //    There is no fallback chain: the lookup is a single lexical hit under the module dir.
    def 'should ignore a same-named skill beside the including script'() {
        given:
        final root = Files.createTempDirectory('test')
        final mod = root.resolve('mods/reporter')
        write(mod.resolve('main.nf'), agentModule('reporter', "skills 'qa-report'"))
        skill(mod.resolve('skills/qa-report'), 'qa-report', MODULE_MARKER)
        and: 'a decoy skill of the SAME name next to the entry script'
        skill(root.resolve('skills/qa-report'), 'qa-report', DECOY_MARKER)
        final main = write(root.resolve('main.nf'),
            entryScript("include { reporter as qc } from './mods/reporter'", 'qc'))

        when:
        runScript(main)

        then:
        captured.skills*.name == ['qa-report']
        captured.skills[0].content.contains(MODULE_MARKER)
        !captured.skills*.content.join('').contains(DECOY_MARKER)
    }

    // -- two DIFFERENT modules may each ship a skill of the same name: no conflict, and each
    //    agent's request carries its OWN module's content (the duplicate-name rejection is
    //    within one agent only)
    def 'should give each module agent its own copy of a same-named skill'() {
        given:
        final root = Files.createTempDirectory('test')
        final modA = root.resolve('mods/alpha')
        write(modA.resolve('main.nf'), agentModule('alpha', "skills 'qa-report'"))
        skill(modA.resolve('skills/qa-report'), 'qa-report', 'ALPHA-BODY')
        final modB = root.resolve('mods/beta')
        write(modB.resolve('main.nf'), agentModule('beta', "skills 'qa-report'"))
        skill(modB.resolve('skills/qa-report'), 'qa-report', 'BETA-BODY')
        final main = write(root.resolve('main.nf'), '''\
            include { alpha } from './mods/alpha'
            include { beta } from './mods/beta'

            workflow {
                alpha(channel.of('s1')).view()
                beta(channel.of('s2')).view()
            }
            '''.stripIndent())

        and: 'capture every request by agent name'
        final Map<String,AgentRunnerRequest> byName = Collections.synchronizedMap([:])
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req ->
            byName.put(req.agentName, req)
            return 'ok'
        } as AgentRunner

        when:
        runScript(main)

        then:
        byName.keySet() == ['alpha', 'beta'] as Set
        byName.alpha.skills*.name == ['qa-report']
        byName.alpha.skills[0].content.contains('ALPHA-BODY')
        !byName.alpha.skills[0].content.contains('BETA-BODY')
        byName.beta.skills[0].content.contains('BETA-BODY')
        !byName.beta.skills[0].content.contains('ALPHA-BODY')
    }

    // -- the failure names the resolved MODULE dir, so the user can see where it looked
    def 'should fail naming the module directory when a declared skill is missing'() {
        given:
        final root = Files.createTempDirectory('test')
        final mod = root.resolve('mods/reporter')
        write(mod.resolve('main.nf'), agentModule('reporter', "skills 'nope'"))
        and: 'a skill of that name exists ONLY beside the entry script'
        skill(root.resolve('skills/nope'), 'nope', DECOY_MARKER)
        final main = write(root.resolve('main.nf'),
            entryScript("include { reporter } from './mods/reporter'", 'reporter'))

        when:
        runScript(main)

        then:
        final e = thrown(Exception)
        final msg = messages(e)
        msg.contains('Agent skill `nope` not found')
        and: 'the reported directory is the MODULE skills dir, not the launch dir'
        msg.contains('mods/reporter/skills/nope')
        !msg.contains("${root.fileName}/skills/nope".toString())
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
