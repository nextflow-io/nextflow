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

package nextflow.lineage.serde

import java.nio.file.Path
import java.time.OffsetDateTime

import nextflow.lineage.model.v1beta1.Checksum
import nextflow.lineage.model.v1beta1.DataPath
import nextflow.lineage.model.v1beta1.FileOutput
import nextflow.lineage.model.v1beta1.Parameter
import nextflow.lineage.model.v1beta1.TaskOutput
import nextflow.lineage.model.v1beta1.TaskRun
import nextflow.lineage.model.v1beta1.Workflow
import nextflow.lineage.model.v1beta1.WorkflowOutput
import nextflow.lineage.model.v1beta1.WorkflowRun
import spock.lang.Specification

class LinEncoderTest extends Specification{

    def 'should encode and decode Outputs'(){
        given:
            def encoder = new LinEncoder()
        and:
            def output = new FileOutput("/path/to/file", new Checksum("hash_value", "hash_algorithm", "standard"),
                "lid://source", "lid://workflow", "lid://task", 1234)

        when:
            def encoded = encoder.encode(output)
            def object = encoder.decode(encoded)

        then:
            object instanceof FileOutput
            def result = object as FileOutput
            result.path == "/path/to/file"
            result.checksum instanceof Checksum
            result.checksum.value == "hash_value"
            result.checksum.algorithm == "hash_algorithm"
            result.checksum.mode == "standard"
            result.source == "lid://source"
            result.size == 1234

    }

    def 'should encode and decode WorkflowRuns'(){
        given:
        def encoder = new LinEncoder()
        and:
        def config = [
            process: [
                container     : "quay.io/nextflow/bash",
                executor      : "local",
                resourceLabels: ["owner": "xxx"],
                scratch       : false
            ]
        ]
        def metadata = [
            runName        : "big_kare",
            start          : "2025-11-06T13:35:42.049135334Z",
            container      : "quay.io/nextflow/bash",
            commandLine    : "nextflow run 'https://github.com/nextflow-io/hello' -name big_kare -with-tower -r master",
            nextflow       : [version: "25.10.0", enable: [dsl: 2.0]],
            containerEngine: "docker",
            wave           : [enabled: true],
            fusion         : [enabled: true, version: "2.4"],
            seqeraPlatform : [
                workflowId: "wf1234",
                user      : [id: "xxx", userName: "john-smith", email: "john.smith@acme.com", firstName: "John", lastName: "Smith", organization: "acme"],
                workspace : [id: "1234", name: "test-workspace", fullName: "Test workspace", organization: "acme"],
                computeEnv: [id: "ce3456", name: "test-ce", platform: "aws-cloud"],
                pipeline  : [id: "pipe294", name: "https://github.com/nextflow-io/hello", revision: "master", commitId: null],
                labels    : []
            ],
            failOnIgnore   : false
        ]
        def uniqueId = UUID.randomUUID()
        def mainScript = new DataPath("file://path/to/main.nf", new Checksum("78910", "nextflow", "standard"))
        def workflow = new Workflow([mainScript], "https://nextflow.io/nf-test/", "123456")
        def wfRun = new WorkflowRun(workflow, uniqueId.toString(), "test_run",
            [new Parameter("String", "param1", "value1"), new Parameter("String", "param2", "value2")],
            config, metadata
        )

        when:
        def encoded = encoder.encode(wfRun)
        def object = encoder.decode(encoded)

        then:
        object instanceof WorkflowRun
        def result = object as WorkflowRun
        result.workflow instanceof Workflow
        result.workflow.scriptFiles.first instanceof DataPath
        result.workflow.scriptFiles.first.path == "file://path/to/main.nf"
        result.workflow.scriptFiles.first.checksum instanceof Checksum
        result.workflow.scriptFiles.first.checksum.value == "78910"
        result.workflow.commitId == "123456"
        result.sessionId == uniqueId.toString()
        result.name == "test_run"
        result.params.size() == 2
        result.params.get(0).name == "param1"
        result.config.process.container == "quay.io/nextflow/bash"
        result.config.process.executor == "local"
        result.metadata.seqeraPlatform.workflowId == "wf1234"

    }

    def 'should decode WorkflowRuns without metadata'(){
        given:
        def encoder = new LinEncoder()
        def wfRunStr = '''
{
  "version": "lineage/v1beta1",
  "kind": "WorkflowRun",
  "spec": {
    "workflow": {
      "scriptFiles": [
        {
          "path": "https://github.com/nextflow-io/hello/main.nf",
          "checksum": {
            "value": "78910",
            "algorithm": "nextflow",
            "mode": "standard"
          }
        }
      ],
      "repository": "https://github.com/nextflow-io/hello",
      "commitId": "2ce0b0e2943449188092a0e25102f0dadc70cb0a"
    },
    "sessionId": "4f02559e-9ebd-41d8-8ee2-a8d1e4f09c67",
    "name": "test_run",
    "params": [],
    "config": {
      "process": {
        "container": "quay.io/nextflow/bash",
        "executor": "local",
        "resourceLabels": {
          "owner": "xxx"
        },
        "scratch": false
      }
    }
  }
}
        '''
        when:
        def object = encoder.decode(wfRunStr)

        then:
        object instanceof WorkflowRun
        def result = object as WorkflowRun
        result.workflow instanceof Workflow
        result.workflow.scriptFiles.first instanceof DataPath
        result.workflow.scriptFiles.first.path == "https://github.com/nextflow-io/hello/main.nf"
        result.workflow.scriptFiles.first.checksum instanceof Checksum
        result.workflow.scriptFiles.first.checksum.value == "78910"
        result.workflow.commitId == "2ce0b0e2943449188092a0e25102f0dadc70cb0a"
        result.sessionId == "4f02559e-9ebd-41d8-8ee2-a8d1e4f09c67"
        result.name == "test_run"
        result.params.size() == 0
        result.config.process.container == "quay.io/nextflow/bash"
        result.config.process.executor == "local"
        result.metadata == null
    }

    def 'should encode and decode WorkflowOutputs'(){
        given:
        def encoder = new LinEncoder()
        and:
        def time = OffsetDateTime.now()
        def wfResults = new WorkflowOutput(time, "lid://1234", [
            new Parameter("Collection", "a", [id: 'id', file: 'sample.txt' as Path]),
            new Parameter("String", "b", "B")
        ])

        when:
        def encoded = encoder.encode(wfResults)
        def object = encoder.decode(encoded)
        then:
        object instanceof WorkflowOutput
        def result = object as WorkflowOutput
        result.createdAt == time
        result.workflowRun == "lid://1234"
        result.output == [new Parameter("Collection", "a", [id: 'id', file: 'sample.txt']), new Parameter("String", "b", "B")]
    }

    def 'should encode and decode TaskRun'() {
        given:
        def encoder = new LinEncoder()
        and:
        def uniqueId = UUID.randomUUID()
        def taskRun = new TaskRun(
            uniqueId.toString(),"name", new Checksum("78910", "nextflow", "standard"), 'this is a script',
            [new Parameter("String", "param1", "value1")], "container:version", "conda", "spack", "amd64",
            [a: "A", b: "B"], [new DataPath("path/to/file", new Checksum("78910", "nextflow", "standard"))]
        )
        when:
        def encoded = encoder.encode(taskRun)
        def object = encoder.decode(encoded)
        then:
        object instanceof TaskRun
        def result = object as TaskRun
        result.sessionId == uniqueId.toString()
        result.name == "name"
        result.codeChecksum.value == "78910"
        result.script == "this is a script"
        result.input.size() == 1
        result.input.get(0).name == "param1"
        result.container == "container:version"
        result.conda == "conda"
        result.spack == "spack"
        result.architecture == "amd64"
        result.globalVars == [a: "A", b: "B"]
        result.binEntries.size() == 1
        result.binEntries.get(0).path == "path/to/file"
        result.binEntries.get(0).checksum.value == "78910"
    }

    def 'should encode and decode TaskOutputs'(){
        given:
        def encoder = new LinEncoder()
        and:
        def time = OffsetDateTime.now()
        def parameter = new Parameter("a","b", "c")
        def wfResults = new TaskOutput("lid://1234", "lid://5678", time, [parameter], null)
        when:
        def encoded = encoder.encode(wfResults)
        def object = encoder.decode(encoded)

        then:
        object instanceof TaskOutput
        def result = object as TaskOutput
        result.createdAt == time
        result.taskRun == "lid://1234"
        result.workflowRun == "lid://5678"
        result.output.size() == 1
        result.output[0] == parameter
    }

    def 'object with null date attributes' () {
        given:
        def encoder = new LinEncoder()
        and:
        def wfResults = new WorkflowOutput(null, "lid://1234")
        when:
        def encoded = encoder.encode(wfResults)
        def object = encoder.decode(encoded)
        then:
        encoded == '{"version":"lineage/v1beta1","kind":"WorkflowOutput","spec":{"createdAt":null,"workflowRun":"lid://1234","output":null}}'
        def result = object as WorkflowOutput
        result.createdAt == null

    }
}
