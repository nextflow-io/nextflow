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
import java.nio.file.attribute.BasicFileAttributes

import nextflow.SysEnv
import nextflow.config.ConfigMap
import nextflow.exception.AbortOperationException
import nextflow.file.FileHelper
import nextflow.lineage.cli.LinCommandImpl
import nextflow.lineage.fs.LinFileSystemProvider
import nextflow.lineage.model.v1beta1.Checksum
import nextflow.lineage.model.v1beta1.DataPath
import nextflow.lineage.model.v1beta1.FileOutput
import nextflow.lineage.model.v1beta1.Parameter
import nextflow.lineage.model.v1beta1.Workflow
import nextflow.lineage.model.v1beta1.WorkflowRun
import nextflow.lineage.serde.LinEncoder
import nextflow.plugin.Plugins
import nextflow.util.CacheHelper
import org.junit.Rule
import spock.lang.Narrative
import spock.lang.See
import spock.lang.Shared
import spock.lang.Specification
import spock.lang.Subject
import spock.lang.Title
import test.OutputCapture

/**
 * Integration tests for lineage validation functionality.
 * Tests end-to-end validation of workflow runs using the CLI command.
 *
 * @author Edmund Miller <edmund.a.miller@gmail.com>
 */
@Title("nextflow lineage validate — end-to-end semantic equivalence check")
@Narrative('''
Drives `LinCommandImpl.validate` against a real on-disk lineage store populated
with WorkflowRun + FileOutput records. These specs exercise the CI-led primary
use case from the ADR: two runs that differ only in ephemeral fields (session,
name, timestamps, absolute paths) must be reported as semantically equivalent;
two runs whose published outputs, parameters, or workflow identity disagree
must abort the command non-zero.

Whenever a behaviour here disagrees with the Spock LineageSnapshotter integration
the shared `LineageValidator` core has drifted — both surfaces must stay in lockstep.
''')
@See([
    "https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md",
    "https://github.com/nextflow-io/nextflow/pull/7167",
    "https://github.com/nextflow-io/nextflow/pull/7168"
])
@Subject(LinCommandImpl)
class LinValidateIntegrationTest extends Specification {

    @Shared
    Path tmpDir

    @Shared
    Path storeLocation

    @Shared
    ConfigMap configMap

    LinEncoder encoder

    def reset() {
        def provider = FileHelper.getProviderFor('lid') as LinFileSystemProvider
        provider?.reset()
        LinStoreFactory.reset()
    }

    def setup() {
        reset()
        SysEnv.push([:])
        tmpDir = Files.createTempDirectory('validate-test')
        storeLocation = tmpDir.resolve("store")
        Files.createDirectories(storeLocation)
        configMap = new ConfigMap([lineage: [enabled: true, store: [location: storeLocation.toString(), logLocation: storeLocation.resolve(".log").toString()]]])
        encoder = new LinEncoder()
    }

    def cleanup() {
        Plugins.stop()
        LinStoreFactory.reset()
        SysEnv.pop()
        tmpDir?.deleteDir()
    }

    def setupSpec() {
        reset()
    }

    @Rule
    OutputCapture capture = new OutputCapture()

    @See([
        "https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md#d2--equivalence-unit-outputs--key-inputs",
        "https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md#d3--file-equivalence-checksum-only"
    ])
    def 'should validate two equivalent workflow runs with outputs'() {
        given: 'Create two workflow runs with same structure'
        def wf = new Workflow([new DataPath("file:///path/to/main.nf")], "https://github.com/test/repo", "abc123")
        
        // First workflow run
        def run1 = new WorkflowRun(wf, "session-1", "eager_einstein",
            [new Parameter("String", "input", "sample.fastq"),
             new Parameter("Integer", "threads", 4)])
        
        // Second workflow run - different session/name but same params
        def run2 = new WorkflowRun(wf, "session-2", "angry_turing",
            [new Parameter("String", "input", "sample.fastq"),
             new Parameter("Integer", "threads", 4)])

        and: 'Create output files with same checksums'
        def outputDir = tmpDir.resolve('outputs')
        Files.createDirectories(outputDir)
        
        def outFile1 = outputDir.resolve('result.txt')
        outFile1.text = 'identical output content'
        def checksum = CacheHelper.hasher(outFile1).hash().toString()
        def attrs = Files.readAttributes(outFile1, BasicFileAttributes)

        and: 'Store workflow runs'
        def lid1 = storeLocation.resolve("wf1/.data.json")
        def lid2 = storeLocation.resolve("wf2/.data.json")
        Files.createDirectories(lid1.parent)
        Files.createDirectories(lid2.parent)
        lid1.text = encoder.encode(run1)
        lid2.text = encoder.encode(run2)

        and: 'Store outputs with same checksums but different timestamps/paths'
        def output1 = new FileOutput(
            outputDir.resolve('run1/result.txt').toString(),
            new Checksum(checksum, "nextflow", "standard"),
            "lid://wf1/result.txt", "lid://wf1", null,
            attrs.size(), LinUtils.toDate(attrs.creationTime()), LinUtils.toDate(attrs.lastModifiedTime())
        )
        def output2 = new FileOutput(
            outputDir.resolve('run2/result.txt').toString(),
            new Checksum(checksum, "nextflow", "standard"),  // Same checksum!
            "lid://wf2/result.txt", "lid://wf2", null,
            attrs.size(), LinUtils.toDate(attrs.creationTime()).plusSeconds(100), // Different time
            LinUtils.toDate(attrs.lastModifiedTime()).plusSeconds(100)
        )
        
        def out1Path = storeLocation.resolve("wf1/result.txt/.data.json")
        def out2Path = storeLocation.resolve("wf2/result.txt/.data.json")
        Files.createDirectories(out1Path.parent)
        Files.createDirectories(out2Path.parent)
        out1Path.text = encoder.encode(output1)
        out2Path.text = encoder.encode(output2)

        when: 'Validate workflow runs'
        new LinCommandImpl().validate(configMap, ["lid://wf1", "--against", "lid://wf2"])
        def stdout = capture
            .toString()
            .readLines()
            .findResults { line -> !line.contains('DEBUG') ? line : null }
            .findResults { line -> !line.contains('INFO') ? line : null }
            .findResults { line -> !line.contains('plugin') ? line : null }

        then: 'Should pass - runs are semantically equivalent'
        noExceptionThrown()
        stdout.any { it.contains("semantically equivalent") }
    }

    @See("https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md#d13--difference-categories-outputs--params--workflow-identity--resources")
    def 'should detect different output checksums between runs'() {
        given: 'Create two workflow runs'
        def wf = new Workflow([new DataPath("file:///path/to/main.nf")], "https://github.com/test/repo", "abc123")
        def run1 = new WorkflowRun(wf, "session-1", "run1", [])
        def run2 = new WorkflowRun(wf, "session-2", "run2", [])

        and: 'Store workflow runs'
        def lid1 = storeLocation.resolve("wf1/.data.json")
        def lid2 = storeLocation.resolve("wf2/.data.json")
        Files.createDirectories(lid1.parent)
        Files.createDirectories(lid2.parent)
        lid1.text = encoder.encode(run1)
        lid2.text = encoder.encode(run2)

        and: 'Store outputs with DIFFERENT checksums'
        def output1 = new FileOutput(
            "/results/out.txt",
            new Checksum("checksum_version_1", "nextflow", "standard"),
            "lid://wf1/out.txt", "lid://wf1", null, 100, null, null
        )
        def output2 = new FileOutput(
            "/results/out.txt",
            new Checksum("checksum_version_2_DIFFERENT", "nextflow", "standard"),  // Different!
            "lid://wf2/out.txt", "lid://wf2", null, 100, null, null
        )
        
        def out1Path = storeLocation.resolve("wf1/out.txt/.data.json")
        def out2Path = storeLocation.resolve("wf2/out.txt/.data.json")
        Files.createDirectories(out1Path.parent)
        Files.createDirectories(out2Path.parent)
        out1Path.text = encoder.encode(output1)
        out2Path.text = encoder.encode(output2)

        when: 'Validate workflow runs'
        new LinCommandImpl().validate(configMap, ["lid://wf1", "--against", "lid://wf2"])

        then: 'Should fail - checksums differ'
        def error = thrown(AbortOperationException)
        error.message.contains("not equivalent")
    }

    @See("https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md#d13--difference-categories-outputs--params--workflow-identity--resources")
    def 'should detect different parameters between runs'() {
        given: 'Create two workflow runs with different params'
        def wf = new Workflow([new DataPath("file:///path/to/main.nf")], "https://github.com/test/repo", "abc123")
        def run1 = new WorkflowRun(wf, "session-1", "run1",
            [new Parameter("String", "input", "sample_A.fastq")])
        def run2 = new WorkflowRun(wf, "session-2", "run2",
            [new Parameter("String", "input", "sample_B.fastq")])  // Different input!

        and: 'Store workflow runs'
        def lid1 = storeLocation.resolve("wf1/.data.json")
        def lid2 = storeLocation.resolve("wf2/.data.json")
        Files.createDirectories(lid1.parent)
        Files.createDirectories(lid2.parent)
        lid1.text = encoder.encode(run1)
        lid2.text = encoder.encode(run2)

        when: 'Validate workflow runs'
        new LinCommandImpl().validate(configMap, ["lid://wf1", "--against", "lid://wf2"])

        then: 'Should fail - params differ'
        def error = thrown(AbortOperationException)
        error.message.contains("not equivalent")
    }

    @See("https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md#d5--subworkflow-handling-flatten-to-terminal-outputs")
    def 'should detect different number of outputs'() {
        given: 'Create two workflow runs'
        def wf = new Workflow([new DataPath("file:///path/to/main.nf")], "https://github.com/test/repo", "abc123")
        def run1 = new WorkflowRun(wf, "session-1", "run1", [])
        def run2 = new WorkflowRun(wf, "session-2", "run2", [])

        and: 'Store workflow runs'
        def lid1 = storeLocation.resolve("wf1/.data.json")
        def lid2 = storeLocation.resolve("wf2/.data.json")
        Files.createDirectories(lid1.parent)
        Files.createDirectories(lid2.parent)
        lid1.text = encoder.encode(run1)
        lid2.text = encoder.encode(run2)

        and: 'First run has 2 outputs, second has 1'
        def output1a = new FileOutput("/results/a.txt", new Checksum("hash_a", "nextflow", "standard"),
            "lid://wf1/a.txt", "lid://wf1", null, 100, null, null)
        def output1b = new FileOutput("/results/b.txt", new Checksum("hash_b", "nextflow", "standard"),
            "lid://wf1/b.txt", "lid://wf1", null, 100, null, null)
        def output2a = new FileOutput("/results/a.txt", new Checksum("hash_a", "nextflow", "standard"),
            "lid://wf2/a.txt", "lid://wf2", null, 100, null, null)
        // Note: wf2 is missing b.txt output
        
        def out1aPath = storeLocation.resolve("wf1/a.txt/.data.json")
        def out1bPath = storeLocation.resolve("wf1/b.txt/.data.json")
        def out2aPath = storeLocation.resolve("wf2/a.txt/.data.json")
        Files.createDirectories(out1aPath.parent)
        Files.createDirectories(out1bPath.parent)
        Files.createDirectories(out2aPath.parent)
        out1aPath.text = encoder.encode(output1a)
        out1bPath.text = encoder.encode(output1b)
        out2aPath.text = encoder.encode(output2a)

        when: 'Validate workflow runs'
        new LinCommandImpl().validate(configMap, ["lid://wf1", "--against", "lid://wf2"])

        then: 'Should fail - different number of outputs'
        def error = thrown(AbortOperationException)
        error.message.contains("not equivalent")
    }

    @See("https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md#d14--output-join-key-relative-path-under--outputdir")
    def 'should validate runs with nested output directories'() {
        given: 'Create workflow runs'
        def wf = new Workflow([new DataPath("file:///path/to/main.nf")], "https://github.com/test/repo", "abc123")
        def run1 = new WorkflowRun(wf, "session-1", "run1", [])
        def run2 = new WorkflowRun(wf, "session-2", "run2", [])

        and: 'Store workflow runs'
        def lid1 = storeLocation.resolve("wf1/.data.json")
        def lid2 = storeLocation.resolve("wf2/.data.json")
        Files.createDirectories(lid1.parent)
        Files.createDirectories(lid2.parent)
        lid1.text = encoder.encode(run1)
        lid2.text = encoder.encode(run2)

        and: 'Create nested outputs with same structure and checksums'
        def paths = ['sample1/aligned.bam', 'sample1/aligned.bam.bai', 'sample2/aligned.bam', 'sample2/aligned.bam.bai']
        paths.each { relativePath ->
            def checksum = "hash_${relativePath.replace('/', '_')}"
            def output1 = new FileOutput("/results/${relativePath}",
                new Checksum(checksum, "nextflow", "standard"),
                "lid://wf1/${relativePath}", "lid://wf1", null, 100, null, null)
            def output2 = new FileOutput("/results/${relativePath}",
                new Checksum(checksum, "nextflow", "standard"),
                "lid://wf2/${relativePath}", "lid://wf2", null, 100, null, null)
            
            def out1Path = storeLocation.resolve("wf1/${relativePath}/.data.json")
            def out2Path = storeLocation.resolve("wf2/${relativePath}/.data.json")
            Files.createDirectories(out1Path.parent)
            Files.createDirectories(out2Path.parent)
            out1Path.text = encoder.encode(output1)
            out2Path.text = encoder.encode(output2)
        }

        when: 'Validate workflow runs'
        new LinCommandImpl().validate(configMap, ["lid://wf1", "--against", "lid://wf2"])
        def stdout = capture
            .toString()
            .readLines()
            .findResults { line -> !line.contains('DEBUG') ? line : null }
            .findResults { line -> !line.contains('INFO') ? line : null }
            .findResults { line -> !line.contains('plugin') ? line : null }

        then: 'Should pass - all nested outputs match'
        noExceptionThrown()
        stdout.any { it.contains("semantically equivalent") }
    }

    @See("https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md#d8--ignore-mechanism-flag--config-jsonpath-style")
    def 'should ignore specified fields during validation'() {
        given: 'Create workflow runs with different commitIds'
        def wf1 = new Workflow([new DataPath("file:///path/to/main.nf")], "https://github.com/test/repo", "commit_v1")
        def wf2 = new Workflow([new DataPath("file:///path/to/main.nf")], "https://github.com/test/repo", "commit_v2")
        def run1 = new WorkflowRun(wf1, "session-1", "run1", [])
        def run2 = new WorkflowRun(wf2, "session-2", "run2", [])

        and: 'Store workflow runs'
        def lid1 = storeLocation.resolve("wf1/.data.json")
        def lid2 = storeLocation.resolve("wf2/.data.json")
        Files.createDirectories(lid1.parent)
        Files.createDirectories(lid2.parent)
        lid1.text = encoder.encode(run1)
        lid2.text = encoder.encode(run2)

        when: 'Validate with --ignore-fields commitId'
        new LinCommandImpl().validate(configMap, 
            ["lid://wf1", "--against", "lid://wf2", "--ignore-fields", "commitId"])
        def stdout = capture
            .toString()
            .readLines()
            .findResults { line -> !line.contains('DEBUG') ? line : null }
            .findResults { line -> !line.contains('INFO') ? line : null }
            .findResults { line -> !line.contains('plugin') ? line : null }

        then: 'Should pass - commitId is ignored'
        noExceptionThrown()
        stdout.any { it.contains("semantically equivalent") }
    }

    @See("https://github.com/nextflow-io/nextflow/blob/master/adr/20260521-lineage-validate.md#d6--failure-output-human-diff-default---json-for-ci")
    def 'should show diff when validation fails'() {
        given: 'Create workflow runs with different configs'
        def wf = new Workflow([new DataPath("file:///path/to/main.nf")], "https://github.com/test/repo", "abc123")
        def run1 = new WorkflowRun(wf, "session-1", "run1",
            [new Parameter("Integer", "memory", 8)])
        def run2 = new WorkflowRun(wf, "session-2", "run2",
            [new Parameter("Integer", "memory", 16)])  // Different memory!

        and: 'Store workflow runs'
        def lid1 = storeLocation.resolve("wf1/.data.json")
        def lid2 = storeLocation.resolve("wf2/.data.json")
        Files.createDirectories(lid1.parent)
        Files.createDirectories(lid2.parent)
        lid1.text = encoder.encode(run1)
        lid2.text = encoder.encode(run2)

        when: 'Validate workflow runs'
        try {
            new LinCommandImpl().validate(configMap, ["lid://wf1", "--against", "lid://wf2"])
        } catch (AbortOperationException e) {
            // Expected
        }
        def stdout = capture.toString()

        then: 'Should show diff output'
        stdout.contains("differ") || stdout.contains("---") || stdout.contains("+++")
    }
}
