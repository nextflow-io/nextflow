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
package nextflow.processor.hash

import java.nio.file.Paths

import nextflow.Session
import nextflow.processor.Architecture
import nextflow.processor.TaskConfig
import nextflow.processor.TaskHasher
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import nextflow.script.ProcessConfig
import nextflow.script.TaskClosure
import nextflow.script.bundle.ResourcesBundle
import nextflow.script.params.InParam
import spock.lang.Specification

class SpecTaskHasherTest extends Specification {

    /**
     * A task exercising every key the unit level can reach.
     *
     * Note: TaskRun.getCondaEnv()/getSpackEnv() return Path and
     * TaskConfig.getArchitecture() returns Architecture, not String — the mock
     * stubs use real instances of those types so Spock's response coercion does
     * not throw a GroovyCastException.
     */
    private Map fixture() {
        def session = Mock(Session) {
            getUniqueId() >> UUID.fromString('b69b6eeb-b332-4d2c-9957-c291b15f498c')
            enableModuleBinaries() >> false
            getStubRun() >> false
        }
        def processor = Mock(TaskProcessor) {
            getName() >> 'PIPE:FOO'
            getSession() >> session
            getConfig() >> Mock(ProcessConfig) { getHashMode() >> nextflow.util.CacheHelper.HashMode.STANDARD }
            getModuleBundle() >> null
        }
        def config = Mock(TaskConfig) {
            getModule() >> ['gcc/11']
            getArchitecture() >> new Architecture('linux/amd64')
            getStubBlock() >> null
        }
        def task = Mock(TaskRun) {
            getSource() >> 'echo hello'
            getProcessor() >> processor
            getConfig() >> config
            isContainerEnabled() >> false
            getInputs() >> [:]
            getOutputEvals() >> [alpha: 'echo a']
            getCondaEnv() >> Paths.get('env.yml')
            getSpackEnv() >> Paths.get('spack.yaml')
        }
        def helper = Spy(new TaskHasher(task))
        helper.getTaskGlobalVars() >> [foo: 'a', bar: 'b']
        helper.getTaskBinEntries(_) >> []
        return [task: task, helper: helper]
    }

    /**
     * Same as {@link #fixture()} but with a stub run enabled: `session.stubRun` is
     * `true`, `task.hasStubBlock()` is `true` (what Contributors.STUB_MARKER reads),
     * and `config.getStubBlock()` returns a non-null TaskClosure (what
     * TaskHasher.compute() reads directly). This covers the STUB_MARKER branch that
     * fixture() leaves dark.
     *
     * Note: under @CompileStatic both TaskHasher.compute() and Contributors.STUB_MARKER
     * resolve `session.stubRun` to `isStubRun()`, not `getStubRun()` — so the mock
     * stubs that method explicitly (see ContributorsTest for the same finding).
     */
    private Map stubRunFixture() {
        def session = Mock(Session) {
            getUniqueId() >> UUID.fromString('b69b6eeb-b332-4d2c-9957-c291b15f498c')
            enableModuleBinaries() >> false
            isStubRun() >> true
        }
        def processor = Mock(TaskProcessor) {
            getName() >> 'PIPE:FOO'
            getSession() >> session
            getConfig() >> Mock(ProcessConfig) { getHashMode() >> nextflow.util.CacheHelper.HashMode.STANDARD }
            getModuleBundle() >> null
        }
        def config = Mock(TaskConfig) {
            getModule() >> ['gcc/11']
            getArchitecture() >> new Architecture('linux/amd64')
            getStubBlock() >> Mock(TaskClosure)
        }
        def task = Mock(TaskRun) {
            getSource() >> 'echo hello'
            getProcessor() >> processor
            getConfig() >> config
            isContainerEnabled() >> false
            getInputs() >> [:]
            getOutputEvals() >> [alpha: 'echo a']
            getCondaEnv() >> Paths.get('env.yml')
            getSpackEnv() >> Paths.get('spack.yaml')
            hasStubBlock() >> true
        }
        def helper = Spy(new TaskHasher(task))
        helper.getTaskGlobalVars() >> [foo: 'a', bar: 'b']
        helper.getTaskBinEntries(_) >> []
        return [task: task, helper: helper]
    }

    /**
     * Same as {@link #fixture()} but exercises the three keys that never fire there:
     * CONTAINER, BIN_ENTRIES (with more than one entry, so a wrong-arity contributor
     * -- one value instead of N -- would be caught), and MODULE_BUNDLE.
     */
    private Map richFixture() {
        def session = Mock(Session) {
            getUniqueId() >> UUID.fromString('b69b6eeb-b332-4d2c-9957-c291b15f498c')
            enableModuleBinaries() >> true
            getStubRun() >> false
        }
        def bundle = Mock(ResourcesBundle) {
            hasEntries() >> true
            fingerprint() >> 'bundle-fingerprint'
        }
        def processor = Mock(TaskProcessor) {
            getName() >> 'PIPE:FOO'
            getSession() >> session
            getConfig() >> Mock(ProcessConfig) { getHashMode() >> nextflow.util.CacheHelper.HashMode.STANDARD }
            getModuleBundle() >> bundle
        }
        def config = Mock(TaskConfig) {
            getModule() >> ['gcc/11']
            getArchitecture() >> new Architecture('linux/amd64')
            getStubBlock() >> null
        }
        def inParam = Mock(InParam) {
            getName() >> 'reads'
        }
        def task = Mock(TaskRun) {
            getSource() >> 'echo hello'
            getProcessor() >> processor
            getConfig() >> config
            isContainerEnabled() >> true
            getContainerFingerprint() >> 'sha256:abcdef'
            getInputs() >> [(inParam): 'value.txt']
            getOutputEvals() >> [alpha: 'echo a']
            getCondaEnv() >> Paths.get('env.yml')
            getSpackEnv() >> Paths.get('spack.yaml')
        }
        def helper = Spy(new TaskHasher(task))
        helper.getTaskGlobalVars() >> [foo: 'a', bar: 'b']
        helper.getTaskBinEntries(_) >> [Paths.get('bin/foo.sh'), Paths.get('bin/bar.sh')]
        return [task: task, helper: helper]
    }

    def 'std/v4 reproduces the legacy TaskHasher byte for byte'() {
        given:
        def f = fixture()
        def legacy = f.helper as TaskHasher
        def ctx = new HashContext(f.task as TaskRun, legacy)

        when:
        def expected = legacy.compute()
        def actual = new SpecTaskHasher(ctx, StdSpecs.STD_V4).compute()

        then:
        actual == expected
    }

    def 'std/v4 reproduces the legacy TaskHasher for a stub-run task'() {
        given:
        def f = stubRunFixture()
        def legacy = f.helper as TaskHasher
        def ctx = new HashContext(f.task as TaskRun, legacy)

        when:
        def expected = legacy.compute()
        def actual = new SpecTaskHasher(ctx, StdSpecs.STD_V4).compute()

        then:
        actual == expected
    }

    def 'std/v4 reproduces the legacy TaskHasher when container, bin entries, and module bundle all fire'() {
        given:
        def f = richFixture()
        def legacy = f.helper as TaskHasher
        def ctx = new HashContext(f.task as TaskRun, legacy)

        when:
        def expected = legacy.compute()
        def actual = new SpecTaskHasher(ctx, StdSpecs.STD_V4).compute()

        then:
        actual == expected
    }

    def 'collectKeys reproduces the legacy key list, element for element'() {
        given:
        def f = fixture()
        def ctx = new HashContext(f.task as TaskRun, f.helper as TaskHasher)

        expect:
        new SpecTaskHasher(ctx, StdSpecs.STD_V4).collectKeys() == [
            UUID.fromString('b69b6eeb-b332-4d2c-9957-c291b15f498c'),
            'PIPE:FOO',
            'echo hello',
            'eval_outputs', [alpha: 'echo a'],
            [foo: 'a', bar: 'b'].entrySet(),
            'gcc/11',
            Paths.get('env.yml'),
            Paths.get('spack.yaml'), new Architecture('linux/amd64')
        ]
    }

    def 'byId resolves known specs and rejects unknown ones'() {
        expect:
        StdSpecs.byId('std/v4').is(StdSpecs.STD_V4)

        when:
        StdSpecs.byId('std/nope')
        then:
        thrown(IllegalArgumentException)
    }
}
