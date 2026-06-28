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
package nextflow.script

import java.nio.file.Files

import nextflow.exception.ScriptRuntimeException
import nextflow.module.ModuleInfo
import nextflow.module.ModuleReference
import spock.lang.Specification
import spock.lang.TempDir

class AgentDefTest extends Specification {

    @TempDir
    File tempDir

    private AgentDef makeAgent(BaseScript script, String name) {
        def prompt = new PromptDef({ -> 'hello' }, 'hello')
        return new AgentDef(script, name, [:], [], [], prompt)
    }

    // Helper: invoke the private static recoverModuleRef via Groovy metaprogramming
    private static ModuleReference recoverModuleRef(java.nio.file.Path moduleDir) {
        def m = AgentDef.getDeclaredMethod('recoverModuleRef', java.nio.file.Path)
        m.accessible = true
        // pass as explicit Object[] so null is not misinterpreted as (Object[]) null (zero-arg)
        return (ModuleReference) m.invoke(null, new Object[]{moduleDir})
    }

    def 'should construct an AgentDef with name and content'() {
        given:
        def script = Mock(BaseScript)

        when:
        def agent = makeAgent(script, 'eval_agent')

        then:
        agent.name == 'eval_agent'
        agent.simpleName == 'eval_agent'
        agent.type == 'agent'
    }

    def 'should fail on run() when the agent does not declare exactly one input'() {
        given:
        def script = Mock(BaseScript)
        def agent = makeAgent(script, 'foo') // no inputs declared

        when:
        agent.run(new Object[0])

        then:
        def e = thrown(ScriptRuntimeException)
        e.message.contains('must declare exactly one input')
    }

    def 'should clone with a new name'() {
        given:
        def script = Mock(BaseScript)
        def agent = makeAgent(script, 'foo')

        when:
        def renamed = agent.cloneWithName('bar')

        then:
        renamed instanceof AgentDef
        renamed.name == 'bar'
        agent.name == 'foo' // original untouched
    }

    def 'should expose the goal directive'() {
        given:
        def directives = [model: 'openai/gpt-5-mini', instruction: 'be careful', goal: 'assemble then QC'] as Map<String,Object>
        def prompt = new PromptDef({ -> 'hi' }, 'hi')
        def agent = new AgentDef(Mock(BaseScript), 'a', directives, [], [], prompt)

        expect:
        agent.goal == 'assemble then QC'
        agent.instruction == 'be careful'
    }

    def 'should return null goal when not declared'() {
        given:
        def agent = new AgentDef(Mock(BaseScript), 'a', [model: 'openai/gpt-5-mini'] as Map<String,Object>, [], [],
            new PromptDef({ -> 'hi' }, 'hi'))
        expect:
        agent.goal == null
    }

    static class TestRec implements nextflow.script.types.Record {}

    def 'should expose the skills directive (single, list, none)'() {
        expect:
        new AgentDef(Mock(BaseScript), 'a', [skills: 'greet'] as Map<String,Object>, [], [], new PromptDef({ -> 'h' }, 'h')).skills == ['greet']
        new AgentDef(Mock(BaseScript), 'a', [skills: ['a', 'b']] as Map<String,Object>, [], [], new PromptDef({ -> 'h' }, 'h')).skills == ['a', 'b']
        new AgentDef(Mock(BaseScript), 'a', [:] as Map<String,Object>, [], [], new PromptDef({ -> 'h' }, 'h')).skills == []
    }

    def 'should reject skills combined with a record (structured) output'() {
        given:
        def inp = new AgentBuilder.AgentInput('q', String)
        def out = new AgentBuilder.AgentOutput('a', TestRec)
        def agent = new AgentDef(Mock(BaseScript), 'a', [skills: 'greet'] as Map<String,Object>, [inp], [out], new PromptDef({ -> 'h' }, 'h'))

        when:
        agent.run(['x'] as Object[])

        then:
        def e = thrown(ScriptRuntimeException)
        e.message.contains('skills')
    }

    // -----------------------------------------------------------------------
    // recoverModuleRef unit tests (offline-safe — no network involved)
    // -----------------------------------------------------------------------

    def 'recoverModuleRef: registry install dir WITH marker returns correct ModuleReference'() {
        given: 'a directory tree matching <base>/modules/<scope>/<name> with .module-info marker'
        def base = tempDir.toPath()
        def moduleDir = base.resolve('modules').resolve('nf-core').resolve('skesa')
        Files.createDirectories(moduleDir)
        moduleDir.resolve(ModuleInfo.MODULE_INFO_FILE).text = 'checksum=abc123'

        when:
        def ref = recoverModuleRef(moduleDir)

        then:
        ref != null
        ref.scope == 'nf-core'
        ref.name == 'skesa'
        ref.fullName == 'nf-core/skesa'
    }

    def 'recoverModuleRef: dir WITHOUT marker returns null (local-file include)'() {
        given: 'same layout but NO .module-info marker file'
        def base = tempDir.toPath()
        def moduleDir = base.resolve('modules').resolve('nf-core').resolve('fastqc')
        Files.createDirectories(moduleDir)
        // intentionally no .module-info

        when:
        def ref = recoverModuleRef(moduleDir)

        then:
        ref == null
    }

    def 'recoverModuleRef: null input returns null'() {
        expect:
        recoverModuleRef(null) == null
    }

    def 'recoverModuleRef: dir with marker but NOT under a "modules" grandparent returns null'() {
        given: 'marker present but parent is named something other than "modules"'
        def base = tempDir.toPath()
        def moduleDir = base.resolve('notmodules').resolve('nf-core').resolve('skesa')
        Files.createDirectories(moduleDir)
        moduleDir.resolve(ModuleInfo.MODULE_INFO_FILE).text = 'checksum=abc'

        when:
        def ref = recoverModuleRef(moduleDir)

        then:
        ref == null
    }

    def 'recoverModuleRef: dir with marker but only one parent level returns null'() {
        given: 'marker present but dir has only one ancestor above it (no scope/name split possible)'
        def base = tempDir.toPath()
        // layout: <base>/modules/skesa  — no scope segment
        def moduleDir = base.resolve('modules').resolve('skesa')
        Files.createDirectories(moduleDir)
        moduleDir.resolve(ModuleInfo.MODULE_INFO_FILE).text = 'checksum=abc'

        when:
        // At this layout: moduleDir.parent.fileName = 'modules', moduleDir.parent.parent.fileName ≠ 'modules'
        // So the check "modulesDir.fileName == 'modules'" fails because the grandparent of
        // moduleDir is the base dir (temp dir), not 'modules'. Returns null.
        def ref = recoverModuleRef(moduleDir)

        then:
        // parent = modules dir (fileName='modules'), parent.parent = base (fileName != 'modules')
        // scope = 'modules', modulesDir = base, base.fileName != 'modules' → null
        ref == null
    }

    def 'recoverModuleRef: multi-segment module name (e.g. nf-core/subworkflows/bam_sort_stats_samtools) returns reference'() {
        given: 'a registry install with a multi-level name stored at scope/first_segment/rest'
        // ModuleStorage stores <base>/modules/<scope>/<name> where name can be 'bam_sort_stats_samtools'
        // (nf-core uses flat single-segment names for modules; multi-segment only for subworkflows,
        //  stored as a single directory name). This test covers an arbitrary valid single-dir name.
        def base = tempDir.toPath()
        def moduleDir = base.resolve('modules').resolve('nf-core').resolve('bwa_mem')
        Files.createDirectories(moduleDir)
        moduleDir.resolve(ModuleInfo.MODULE_INFO_FILE).text = 'checksum=xyz'

        when:
        def ref = recoverModuleRef(moduleDir)

        then:
        ref != null
        ref.scope == 'nf-core'
        ref.name == 'bwa_mem'
    }
}
