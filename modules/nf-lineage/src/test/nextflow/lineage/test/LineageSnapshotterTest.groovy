/*
 * Copyright 2013-2025, Seqera Labs
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

package nextflow.lineage.test

import java.nio.file.Files
import java.nio.file.Path
import java.time.Instant
import java.time.OffsetDateTime
import java.time.ZoneOffset

import nextflow.lineage.LinStore
import nextflow.lineage.model.v1beta1.Checksum
import nextflow.lineage.model.v1beta1.DataPath
import nextflow.lineage.model.v1beta1.FileOutput
import nextflow.lineage.model.v1beta1.Parameter
import nextflow.lineage.model.v1beta1.Workflow
import nextflow.lineage.model.v1beta1.WorkflowRun
import nextflow.lineage.serde.LinSerializable
import spock.lang.Narrative
import spock.lang.See
import spock.lang.Specification
import spock.lang.Subject
import spock.lang.TempDir
import spock.lang.Title

/**
 * Tests for LineageSnapshotter
 *
 * @author Edmund Miller <edmund.a.miller@gmail.com>
 */
@Title("Spock snapshot integration for lineage validation")
@Narrative('''
LineageSnapshotter is the Spock-facing wrapper around the same validator core
the CLI uses. It supports two flavours of assertion:

 * `assertMatchesSnapshot(lid, id)` — compare a workflow run against a baseline
   snapshot file on disk. The file is created on first run; subsequent runs
   compare against it. `UPDATE_SNAPSHOTS=true` refreshes the baseline.
 * `assertEquivalent(lidA, lidB)` — compare two runs directly, no file.

Snapshot file format and ignore semantics must stay aligned with the CLI: if a
test passes locally but `nextflow lineage validate` fails in CI (or vice versa)
the abstraction has leaked.
''')
@See([
    "https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md#d4--baseline-model-peer-lid-and-snapshot-file-from-day-one",
    "https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md#d11--cli-vs-lineagesnapshotter-extract-shared-lineagevalidator",
    "https://spockframework.org/spock/docs/2.4/all_in_one.html#_specifications_as_documentation"
])
@Subject(LineageSnapshotter)
class LineageSnapshotterTest extends Specification {

    @TempDir
    Path tmpDir

    Path snapshotDir
    LinStore mockStore

    def setup() {
        snapshotDir = tmpDir.resolve('snapshots')
        Files.createDirectories(snapshotDir)
        mockStore = Mock(LinStore)
    }

    @See("https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md#d4--baseline-model-peer-lid-and-snapshot-file-from-day-one")
    def 'should create snapshot on first run'() {
        given:
        def wf = new Workflow([new DataPath("/path/to/main.nf")], "test-wf", "abc123")
        def run = new WorkflowRun(wf, "session-1", "run1", 
            [new Parameter("String", "input", "test.fastq")])
        
        mockStore.load("wf1") >> run
        mockStore.getSubKeys("wf1") >> { [].stream() }
        
        def snapshotter = new LineageSnapshotter(mockStore)
            .withSnapshotDir(snapshotDir)

        when:
        snapshotter.assertMatchesSnapshot("lid://wf1", "test-baseline")

        then:
        noExceptionThrown()
        Files.exists(snapshotDir.resolve("test-baseline.json"))
    }

    @See("https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md#d2--equivalence-unit-outputs--key-inputs")
    def 'should pass when snapshot matches'() {
        given:
        def wf = new Workflow([new DataPath("/path/to/main.nf")], "test-wf", "abc123")
        def run1 = new WorkflowRun(wf, "session-1", "run1", 
            [new Parameter("String", "input", "test.fastq")])
        def run2 = new WorkflowRun(wf, "session-2", "run2",
            [new Parameter("String", "input", "test.fastq")])
        
        // First call returns run1 (for creating snapshot)
        // Second call returns run2 (for comparing)
        mockStore.load("wf1") >>> [run1, run2]
        mockStore.getSubKeys("wf1") >> { [].stream() }
        
        def snapshotter = new LineageSnapshotter(mockStore)
            .withSnapshotDir(snapshotDir)

        when:
        // First run - creates snapshot
        snapshotter.assertMatchesSnapshot("lid://wf1", "test-baseline")
        // Second run - compares against snapshot
        snapshotter.assertMatchesSnapshot("lid://wf1", "test-baseline")

        then:
        noExceptionThrown()
    }

    @See("https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md#d6--failure-output-human-diff-default---json-for-ci")
    def 'should fail when snapshot differs'() {
        given:
        def wf = new Workflow([new DataPath("/path/to/main.nf")], "test-wf", "abc123")
        def run1 = new WorkflowRun(wf, "session-1", "run1",
            [new Parameter("String", "input", "test.fastq")])
        def run2 = new WorkflowRun(wf, "session-2", "run2",
            [new Parameter("String", "input", "different.fastq")])  // Different param
        
        mockStore.load("wf1") >>> [run1, run2]
        mockStore.getSubKeys("wf1") >> { [].stream() }
        
        def snapshotter = new LineageSnapshotter(mockStore)
            .withSnapshotDir(snapshotDir)

        when:
        // First run - creates snapshot
        snapshotter.assertMatchesSnapshot("lid://wf1", "test-baseline")
        // Second run - should fail
        snapshotter.assertMatchesSnapshot("lid://wf1", "test-baseline")

        then:
        def error = thrown(AssertionError)
        error.message.contains("does not match snapshot")
        error.message.contains("test-baseline")
        // Should create .actual.json for debugging
        Files.exists(snapshotDir.resolve("test-baseline.actual.json"))
    }

    @See("https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md#d4--baseline-model-peer-lid-and-snapshot-file-from-day-one")
    def 'should update snapshot when updateSnapshots is true'() {
        given:
        def wf = new Workflow([new DataPath("/path/to/main.nf")], "test-wf", "abc123")
        def run1 = new WorkflowRun(wf, "session-1", "run1",
            [new Parameter("String", "input", "old.fastq")])
        def run2 = new WorkflowRun(wf, "session-2", "run2",
            [new Parameter("String", "input", "new.fastq")])
        
        mockStore.load("wf1") >>> [run1, run2]
        mockStore.getSubKeys("wf1") >> { [].stream() }
        
        def snapshotter = new LineageSnapshotter(mockStore)
            .withSnapshotDir(snapshotDir)

        when:
        // First run - creates snapshot
        snapshotter.assertMatchesSnapshot("lid://wf1", "test-baseline")
        
        // Enable update mode
        snapshotter.withUpdateSnapshots(true)
        
        // Second run - updates snapshot instead of failing
        snapshotter.assertMatchesSnapshot("lid://wf1", "test-baseline")

        then:
        noExceptionThrown()
        // Snapshot should contain new value
        snapshotDir.resolve("test-baseline.json").text.contains("new.fastq")
    }

    @See("https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md#d11--cli-vs-lineagesnapshotter-extract-shared-lineagevalidator")
    def 'should compare two runs directly with assertEquivalent'() {
        given:
        def wf = new Workflow([new DataPath("/path/to/main.nf")], "test-wf", "abc123")
        def run1 = new WorkflowRun(wf, "session-1", "run1",
            [new Parameter("String", "input", "test.fastq")])
        def run2 = new WorkflowRun(wf, "session-2", "run2",
            [new Parameter("String", "input", "test.fastq")])
        
        mockStore.load("wf1") >> run1
        mockStore.load("wf2") >> run2
        mockStore.getSubKeys(_) >> { [].stream() }
        
        def snapshotter = new LineageSnapshotter(mockStore)
            .withSnapshotDir(snapshotDir)

        when:
        snapshotter.assertEquivalent("lid://wf1", "lid://wf2")

        then:
        noExceptionThrown()
    }

    def 'should fail assertEquivalent when runs differ'() {
        given:
        def wf = new Workflow([new DataPath("/path/to/main.nf")], "test-wf", "abc123")
        def run1 = new WorkflowRun(wf, "session-1", "run1",
            [new Parameter("String", "input", "test.fastq")])
        def run2 = new WorkflowRun(wf, "session-2", "run2",
            [new Parameter("String", "input", "different.fastq")])
        
        mockStore.load("wf1") >> run1
        mockStore.load("wf2") >> run2
        mockStore.getSubKeys(_) >> { [].stream() }
        
        def snapshotter = new LineageSnapshotter(mockStore)
            .withSnapshotDir(snapshotDir)

        when:
        snapshotter.assertEquivalent("lid://wf1", "lid://wf2")

        then:
        def error = thrown(AssertionError)
        error.message.contains("not semantically equivalent")
    }

    def 'should throw if snapshotDir not configured'() {
        given:
        def wf = new Workflow([new DataPath("/path/to/main.nf")], "test-wf", "abc123")
        def run = new WorkflowRun(wf, "session-1", "run1", [])
        
        mockStore.load("wf1") >> run
        mockStore.getSubKeys("wf1") >> { [].stream() }
        
        def snapshotter = new LineageSnapshotter(mockStore)
        // Not calling withSnapshotDir()

        when:
        snapshotter.assertMatchesSnapshot("lid://wf1", "test-baseline")

        then:
        thrown(IllegalStateException)
    }

    def 'should include outputs in snapshot comparison'() {
        given:
        def time = OffsetDateTime.ofInstant(Instant.ofEpochMilli(123456789), ZoneOffset.UTC)
        def wf = new Workflow([new DataPath("/path/to/main.nf")], "test-wf", "abc123")
        def run = new WorkflowRun(wf, "session-1", "run1", [])
        def output = new FileOutput("/results/out.txt", 
            new Checksum("hash123", "nextflow", "standard"),
            "lid://wf1/out.txt", "lid://wf1", null, 1024, time, time, null)
        
        mockStore.load("wf1") >> run
        mockStore.load("wf1/out.txt") >> output
        mockStore.getSubKeys("wf1") >> ["wf1/out.txt"].stream()
        
        def snapshotter = new LineageSnapshotter(mockStore)
            .withSnapshotDir(snapshotDir)

        when:
        snapshotter.assertMatchesSnapshot("lid://wf1", "with-outputs")

        then:
        noExceptionThrown()
        def snapshotContent = snapshotDir.resolve("with-outputs.json").text
        snapshotContent.contains("hash123")
        snapshotContent.contains("out.txt")
    }

    def 'should return normalized JSON for debugging'() {
        given:
        def wf = new Workflow([new DataPath("/path/to/main.nf")], "test-wf", "abc123")
        def run = new WorkflowRun(wf, "session-1", "run1",
            [new Parameter("String", "input", "test.fastq")])
        
        mockStore.load("wf1") >> run
        mockStore.getSubKeys("wf1") >> { [].stream() }
        
        def snapshotter = new LineageSnapshotter(mockStore)
            .withSnapshotDir(snapshotDir)

        when:
        def json = snapshotter.getNormalizedJson("lid://wf1")

        then:
        json.contains("WorkflowRun")
        json.contains("test.fastq")
        !json.contains("session-1")  // sessionId should be stripped
    }

    @See("https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md#d8--ignore-mechanism-flag--config-jsonpath-style")
    def 'should support custom ignore fields'() {
        given:
        def wf1 = new Workflow([new DataPath("/path/to/main.nf")], "test-wf", "commit1")
        def wf2 = new Workflow([new DataPath("/path/to/main.nf")], "test-wf", "commit2")
        def run1 = new WorkflowRun(wf1, "session-1", "run1", [])
        def run2 = new WorkflowRun(wf2, "session-2", "run2", [])
        
        mockStore.load("wf1") >> run1
        mockStore.load("wf2") >> run2
        mockStore.getSubKeys(_) >> { [].stream() }
        
        def snapshotter = new LineageSnapshotter(mockStore)
            .withSnapshotDir(snapshotDir)
            .withIgnoreFields(['commitId'])

        when:
        snapshotter.assertEquivalent("lid://wf1", "lid://wf2")

        then:
        noExceptionThrown()
    }
}
