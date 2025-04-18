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

package nextflow.lineage.serde

import nextflow.lineage.model.Checksum
import nextflow.lineage.model.DataPath
import nextflow.lineage.model.Parameter
import nextflow.lineage.model.DataOutput
import nextflow.lineage.model.TaskOutputs
import nextflow.lineage.model.TaskRun
import nextflow.lineage.model.Workflow
import nextflow.lineage.model.WorkflowOutputs
import nextflow.lineage.model.WorkflowRun
import spock.lang.Specification

import java.time.OffsetDateTime

class LinEncoderTest extends Specification{

    def 'should encode and decode Outputs'(){
        given:
            def encoder = new LinEncoder()
        and:
            def output = new DataOutput("/path/to/file", new Checksum("hash_value", "hash_algorithm", "standard"),
                "lid://source", "lid://workflow", "lid://task", 1234)

        when:
            def encoded = encoder.encode(output)
            def object = encoder.decode(encoded)

        then:
            object instanceof DataOutput
            def result = object as DataOutput
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
        def uniqueId = UUID.randomUUID()
        def mainScript = new DataPath("file://path/to/main.nf", new Checksum("78910", "nextflow", "standard"))
        def workflow = new Workflow([mainScript], "https://nextflow.io/nf-test/", "123456")
        def wfRun = new WorkflowRun(workflow, uniqueId.toString(), "test_run", [new Parameter("String", "param1", "value1"), new Parameter("String", "param2", "value2")])

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
    }

    def 'should encode and decode WorkflowResults'(){
        given:
        def encoder = new LinEncoder()
        and:
        def time = OffsetDateTime.now()
        def wfResults = new WorkflowOutputs(time, "lid://1234", [new Parameter("String", "a", "A"), new Parameter("String", "b", "B")])
        when:
        def encoded = encoder.encode(wfResults)
        def object = encoder.decode(encoded)

        then:
        object instanceof WorkflowOutputs
        def result = object as WorkflowOutputs
        result.createdAt == time
        result.workflowRun == "lid://1234"
        result.outputs == [new Parameter("String", "a", "A"), new Parameter("String", "b", "B")]
    }

    def 'should encode and decode TaskRun'() {
        given:
        def encoder = new LinEncoder()
        and:
        def uniqueId = UUID.randomUUID()
        def taskRun = new TaskRun(
            uniqueId.toString(),"name", new Checksum("78910", "nextflow", "standard"), new Checksum("74517", "nextflow", "standard"),
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
        result.scriptChecksum.value == "74517"
        result.inputs.size() == 1
        result.inputs.get(0).name == "param1"
        result.container == "container:version"
        result.conda == "conda"
        result.spack == "spack"
        result.architecture == "amd64"
        result.globalVars == [a: "A", b: "B"]
        result.binEntries.size() == 1
        result.binEntries.get(0).path == "path/to/file"
        result.binEntries.get(0).checksum.value == "78910"
    }

    def 'should encode and decode TaskResults'(){
        given:
        def encoder = new LinEncoder()
        and:
        def time = OffsetDateTime.now()
        def parameter = new Parameter("a","b", "c")
        def wfResults = new TaskOutputs("lid://1234", "lid://5678", time, [parameter], null)
        when:
        def encoded = encoder.encode(wfResults)
        def object = encoder.decode(encoded)

        then:
        object instanceof TaskOutputs
        def result = object as TaskOutputs
        result.createdAt == time
        result.taskRun == "lid://1234"
        result.workflowRun == "lid://5678"
        result.outputs.size() == 1
        result.outputs[0] == parameter
    }

    def 'object with null date attributes' () {
        given:
        def encoder = new LinEncoder()
        and:
        def wfResults = new WorkflowOutputs(null, "lid://1234")
        when:
        def encoded = encoder.encode(wfResults)
        def object = encoder.decode(encoded)
        then:
        encoded == '{"type":"WorkflowOutputs","createdAt":null,"workflowRun":"lid://1234","outputs":null,"annotations":null}'
        def result = object as WorkflowOutputs
        result.createdAt == null

    }
}
