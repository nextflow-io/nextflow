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

import nextflow.Global
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
import nextflow.util.CacheHelper
import spock.lang.Specification
import spock.lang.Unroll

/**
 * Guards the two invariants that keep published hash versions honest.
 *
 * 1. The newest spec must agree with the default hashing path. They drift apart
 *    the moment someone changes what master hashes, which is exactly when a new
 *    spec is due.
 * 2. Every published spec must keep producing the hash it produced when it was
 *    published. The expected values below are frozen: a change to a contributor,
 *    to {@code HashBuilder} or to the spec files moves them, and an old cache
 *    stops resolving.
 *
 * Both run on one task that fires all thirteen {@link HashKey}s at once, so a
 * change to any key is caught.
 */
class StdSpecsGoldenTest extends Specification {

    /**
     * Frozen hash of {@link #fullTask} under each published spec.
     *
     * Never edit an entry to make the build pass. A moved value means either a
     * deliberate new behaviour — which needs a new spec id, not a new number here —
     * or an accidental change to an old one, which is a cache-breaking bug.
     */
    static final Map<String,String> EXPECTED = [
        'std/v1.1': 'e9cd7f245e4cdaa0db3b5c467c69a693',
        'std/v1.2': 'e9cd7f245e4cdaa0db3b5c467c69a693',
        'std/v1.3': '88107f5e92ffa46b5c275dc26cdb6682',
        'std/v1.4': '59c7644cc45389318476f8b8efc6d833',
        'std/v1.5': 'ac28199a80a5ab576d0d159054a1e59c',
        'std/v1.6': '01b90131cb2a9fe4ea8460ad03938fdd',
        'std/v1.7': 'fdc3e5185e5a69be0c99442730e0c34f',
    ]

    /**
     * Frozen fingerprint of each published spec — a hash of its id, its ordered keys with
     * their contributor names, its encoding rules and its hash function.
     *
     * This catches a changed definition even when the task hash cannot see it: {@code std/v1.1}
     * and {@code std/v1.2} differ only in asset detection, which needs a file inside a Git
     * repository and so produces the same hash in a unit fixture.
     */
    static final Map<String,String> EXPECTED_FINGERPRINT = [
        'std/v1.1': '4af1b4bfdc2157aa95d817068af6c520',
        'std/v1.2': '8f241dc7acc11fd81879dbc252f73b84',
        'std/v1.3': '50bd45e81b4298655d904551cd560fd2',
        'std/v1.4': '73e826e78755090c4adbb64adb3f0a9c',
        'std/v1.5': 'ee3c4b658453611f640f8e6eda403d25',
        'std/v1.6': '51d8389f196c44fdaeb70cc05d9a7ac9',
        'std/v1.7': '90c55a194828689835ec580d12e73057',
    ]

    def setup() {
        // HashBuilder.isAssetFile() reads the global session; leaving another test's
        // session in place would make the file keys depend on test order
        Global.session = null
    }

    def cleanup() {
        Global.session = null
    }

    /**
     * A task firing every key: container, inputs, eval outputs, script vars, bin
     * entries, resources bundle, env modules, conda, spack with architecture, and the
     * stub marker.
     *
     * Every path is absolute. A relative path is resolved against the working
     * directory before hashing, which would make the expected values depend on where
     * the repository is checked out.
     */
    private Map fullTask() {
        def session = Mock(Session) {
            getUniqueId() >> UUID.fromString('b69b6eeb-b332-4d2c-9957-c291b15f498c')
            enableModuleBinaries() >> true
            isStubRun() >> true
        }
        def bundle = Mock(ResourcesBundle) {
            // the contributor tests the bundle for Groovy truth before calling
            // hasEntries(), and an unstubbed asBoolean() on a mock is false
            asBoolean() >> true
            hasEntries() >> true
            fingerprint() >> 'bundle-fingerprint'
        }
        def processor = Mock(TaskProcessor) {
            getName() >> 'PIPE:FOO'
            getSession() >> session
            getConfig() >> Mock(ProcessConfig) { getHashMode() >> CacheHelper.HashMode.STANDARD }
            getModuleBundle() >> bundle
        }
        def config = Mock(TaskConfig) {
            getModule() >> ['gcc/11', 'samtools/1.19']
            getArchitecture() >> new Architecture('linux/amd64')
            getStubBlock() >> Mock(TaskClosure)
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
            // a Map value, so the orderIndependentMaps encoding rule is observable
            getInputs() >> [(inParam): [beta: 2, alpha: 1]]
            getOutputEvals() >> [alpha: 'echo a', beta: 'echo b']
            getCondaEnv() >> Paths.get('/pipeline/env.yml')
            getSpackEnv() >> Paths.get('/pipeline/spack.yaml')
            hasStubBlock() >> true
        }
        def helper = Spy(new TaskHasher(task))
        helper.getTaskGlobalVars() >> [foo: 'a', 'params.bar': 'b']
        helper.getTaskBinEntries(_) >> [Paths.get('/pipeline/bin/foo.sh'), Paths.get('/pipeline/bin/bar.sh')]
        return [task: task, helper: helper]
    }

    private HashContext context(Map f) {
        return new HashContext(f.task as TaskRun, f.helper as TaskHasher)
    }

    def 'the latest spec reproduces the default hashing path'() {
        given:
        def f = fullTask()
        def latest = StdSpecs.latest()

        when:
        def expected = (f.helper as TaskHasher).compute()
        def actual = new VersionedTaskHasher(context(f), latest).compute()

        then:
        assert actual == expected, """\
The default hashing path no longer agrees with ${latest.id}.
Master's task hash changed, so ${latest.id} no longer describes it. Publish a new
spec for the new behaviour instead of editing ${latest.id} -- every cache written
under ${latest.id} still expects the old hash."""
    }

    def 'every key of the latest spec reproduces the default hashing path'() {
        given:
        def f = fullTask()

        expect: 'no key is silently empty, otherwise the test above proves nothing'
        new VersionedTaskHasher(context(f), StdSpecs.latest()).explain().keySet() == HashKey.values() as Set
    }

    @Unroll
    def 'spec #id still produces its published hash'() {
        given:
        def f = fullTask()

        expect:
        assert new VersionedTaskHasher(context(f), StdSpecs.byId(id)).compute().toString() == EXPECTED[id], """\
Spec ${id} no longer produces the hash it was published with.
Something changed how an already published version hashes -- a contributor, HashBuilder,
or the spec file itself. Every cache written under ${id} stops resolving. Find the change
rather than updating the expected value."""

        where:
        id << StdSpecs.IDS
    }

    @Unroll
    def 'spec #id still has its published fingerprint'() {
        expect:
        assert StdSpecs.byId(id).fingerprint().toString() == EXPECTED_FINGERPRINT[id], """\
The definition of spec ${id} changed: ${StdSpecs.byId(id).canonicalForm()}"""

        where:
        id << StdSpecs.IDS
    }

    def 'every published spec is covered by both tables'() {
        expect:
        EXPECTED.keySet() == StdSpecs.IDS as Set
        EXPECTED_FINGERPRINT.keySet() == StdSpecs.IDS as Set

        and: 'no two specs share a fingerprint, otherwise one of them is a duplicate'
        EXPECTED_FINGERPRINT.values().toSet().size() == StdSpecs.IDS.size()
    }
}
