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

package nextflow.lineage

import java.nio.file.Files
import java.nio.file.Path
import java.time.Instant
import java.time.OffsetDateTime
import java.time.ZoneOffset

import nextflow.lineage.model.v1beta1.Checksum
import nextflow.lineage.model.v1beta1.DataPath
import nextflow.lineage.model.v1beta1.FileOutput
import nextflow.lineage.model.v1beta1.Parameter
import nextflow.lineage.model.v1beta1.TaskRun
import nextflow.lineage.model.v1beta1.Workflow
import nextflow.lineage.model.v1beta1.WorkflowOutput
import nextflow.lineage.model.v1beta1.WorkflowRun
import spock.lang.Specification
import spock.lang.TempDir

/**
 * Tests for LinNormalizer
 *
 * @author Edmund Miller <edmund.a.miller@gmail.com>
 */
class LinNormalizerTest extends Specification {

    @TempDir
    Path tmpDir

    def 'should strip ephemeral fields from WorkflowRun'() {
        given:
        def normalizer = new LinNormalizer()
        def wf = new Workflow([new DataPath("/path/to/main.nf")], "hello-nf", "abc123")
        def run = new WorkflowRun(wf, "session-123", "crazy_einstein",
            [new Parameter("String", "input", "test.fastq")])

        when:
        def normalized = normalizer.normalize(run)

        then:
        normalized['kind'] == 'WorkflowRun'
        normalized['spec']['workflow'] != null
        // sessionId should be stripped
        !normalized['spec'].containsKey('sessionId')
        // name should be stripped (auto-generated run name)
        !normalized['spec'].containsKey('name')
        // params should be preserved
        normalized['spec']['params'] != null
    }

    def 'should strip timestamps from FileOutput'() {
        given:
        def normalizer = new LinNormalizer()
        def time = OffsetDateTime.ofInstant(Instant.ofEpochMilli(123456789), ZoneOffset.UTC)
        def output = new FileOutput(
            "/results/sample1/file.bam",
            new Checksum("abc123", "nextflow", "standard"),
            "lid://task123/file.bam",
            "lid://workflow123",
            "lid://task123",
            1024,
            time,
            time,
            ["experiment=test"]
        )

        when:
        def normalized = normalizer.normalize(output)

        then:
        normalized['kind'] == 'FileOutput'
        // createdAt and modifiedAt should be stripped
        !normalized['spec'].containsKey('createdAt')
        !normalized['spec'].containsKey('modifiedAt')
        // checksum should be preserved
        normalized['spec']['checksum']['value'] == 'abc123'
        // labels should be preserved
        normalized['spec']['labels'] == ['experiment=test']
    }

    def 'should strip path field from FileOutput'() {
        given:
        def normalizer = new LinNormalizer()
        def output = new FileOutput(
            "/home/user/results/sample1/file.bam",
            new Checksum("abc123", "nextflow", "standard"),
            "lid://task123/file.bam",
            "lid://workflow123",
            null,
            1024,
            null,
            null,
            null
        )

        when:
        def normalized = normalizer.normalize(output)

        then:
        // path is stripped because absolute paths vary between runs
        // checksum is what matters for comparing outputs
        !normalized['spec'].containsKey('path')
        normalized['spec']['checksum']['value'] == 'abc123'
    }

    def 'should strip path regardless of outputBase setting'() {
        given:
        def normalizer = new LinNormalizer()
            .withOutputBase('/home/user/results')
        def output = new FileOutput(
            "/other/path/file.bam",
            new Checksum("abc123", "nextflow", "standard"),
            "lid://task123/file.bam",
            "lid://workflow123",
            null,
            1024,
            null,
            null,
            null
        )

        when:
        def normalized = normalizer.normalize(output)

        then:
        // path is always stripped as an ephemeral field
        !normalized['spec'].containsKey('path')
    }

    def 'should strip sessionId from TaskRun'() {
        given:
        def normalizer = new LinNormalizer()
        def task = new TaskRun(
            "session-456",
            "ALIGN",
            new Checksum("hash123", "nextflow", "standard"),
            'bwa mem ref.fa reads.fq > out.bam',
            [new Parameter("path", "reads", "reads.fq")],
            "docker://biocontainers/bwa",
            null,
            null,
            null,
            [:],
            []
        )

        when:
        def normalized = normalizer.normalize(task)

        then:
        normalized['kind'] == 'TaskRun'
        !normalized['spec'].containsKey('sessionId')
        normalized['spec']['name'] == null  // name is also stripped
        normalized['spec']['script'] == 'bwa mem ref.fa reads.fq > out.bam'
        normalized['spec']['codeChecksum']['value'] == 'hash123'
    }

    def 'should support additional ignore fields'() {
        given:
        def normalizer = new LinNormalizer()
            .withIgnoreFields(['container', 'conda'])
        def task = new TaskRun(
            "session-456",
            "ALIGN",
            new Checksum("hash123", "nextflow", "standard"),
            'script',
            null,
            "docker://bwa",
            "bioconda::bwa",
            null,
            null,
            [:],
            []
        )

        when:
        def normalized = normalizer.normalize(task)

        then:
        !normalized['spec'].containsKey('container')
        !normalized['spec'].containsKey('conda')
    }

    def 'should compare two identical normalized trees as equal'() {
        given:
        def tree1 = [
            workflowRun: [kind: 'WorkflowRun', spec: [workflow: [name: 'test']]],
            outputs: [
                'file.bam': [kind: 'FileOutput', spec: [checksum: [value: 'abc']]]
            ]
        ]
        def tree2 = [
            workflowRun: [kind: 'WorkflowRun', spec: [workflow: [name: 'test']]],
            outputs: [
                'file.bam': [kind: 'FileOutput', spec: [checksum: [value: 'abc']]]
            ]
        ]

        when:
        def diff = LinNormalizer.compare(tree1, tree2)

        then:
        diff.isEmpty()
    }

    def 'should detect differences in normalized trees'() {
        given:
        def tree1 = [
            workflowRun: [kind: 'WorkflowRun', spec: [workflow: [name: 'test']]],
            outputs: [
                'file.bam': [kind: 'FileOutput', spec: [checksum: [value: 'abc']]]
            ]
        ]
        def tree2 = [
            workflowRun: [kind: 'WorkflowRun', spec: [workflow: [name: 'test']]],
            outputs: [
                'file.bam': [kind: 'FileOutput', spec: [checksum: [value: 'xyz']]]
            ]
        ]

        when:
        def diff = LinNormalizer.compare(tree1, tree2)

        then:
        !diff.isEmpty()
        diff['outputs.file.bam.spec.checksum.value']['expected'] == 'abc'
        diff['outputs.file.bam.spec.checksum.value']['actual'] == 'xyz'
    }

    def 'should detect missing keys'() {
        given:
        def tree1 = [outputs: ['a.txt': [checksum: 'abc'], 'b.txt': [checksum: 'def']]]
        def tree2 = [outputs: ['a.txt': [checksum: 'abc']]]

        when:
        def diff = LinNormalizer.compare(tree1, tree2)

        then:
        !diff.isEmpty()
        diff.containsKey('outputs.b.txt')
    }

    def 'should normalize workflow tree from store'() {
        given:
        def store = Mock(LinStore)
        def normalizer = new LinNormalizer()
        def wf = new Workflow([new DataPath("/path/to/main.nf")], "test-wf", "commit123")
        def run = new WorkflowRun(wf, "session-1", "run_name", [])
        def output = new FileOutput("/results/out.txt", new Checksum("hash", "nf", "std"),
            "lid://wf1/out.txt", "lid://wf1", null, 100, null, null, null)

        store.load("wf1") >> run
        store.getSubKeys("wf1") >> ["wf1/out.txt"].stream()
        store.load("wf1/out.txt") >> output

        when:
        def tree = normalizer.normalizeWorkflowTree(store, "lid://wf1")

        then:
        tree['workflowRun'] != null
        tree['workflowRun']['kind'] == 'WorkflowRun'
        tree['outputs'] != null
        tree['outputs']['out.txt'] != null
        tree['outputs']['out.txt']['kind'] == 'FileOutput'
    }

    def 'should throw when workflow not found'() {
        given:
        def store = Mock(LinStore)
        def normalizer = new LinNormalizer()
        store.load("nonexistent") >> null

        when:
        normalizer.normalizeWorkflowTree(store, "lid://nonexistent")

        then:
        thrown(IllegalArgumentException)
    }
}
