package nextflow.data.cid.serde

import nextflow.data.cid.model.Checksum
import nextflow.data.cid.model.DataPath
import nextflow.data.cid.model.Parameter
import nextflow.data.cid.model.DataOutput
import nextflow.data.cid.model.TaskOutputs
import nextflow.data.cid.model.TaskRun
import nextflow.data.cid.model.Workflow
import nextflow.data.cid.model.WorkflowOutputs
import nextflow.data.cid.model.WorkflowRun
import spock.lang.Specification

import java.time.Instant

class CidEncoderTest extends Specification{

    def 'should encode and decode Outputs'(){
        given:
            def encoder = new CidEncoder()
        and:
            def output = new DataOutput("/path/to/file", new Checksum("hash_value", "hash_algorithm", "standard"), "cid://source", "cid://run", 1234)

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
            result.source == "cid://source"
            result.size == 1234

    }

    def 'should encode and decode WorkflowRuns'(){
        given:
        def encoder = new CidEncoder()
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
        def encoder = new CidEncoder()
        and:
        def time = Instant.now()
        def wfResults = new WorkflowOutputs(time, "cid://1234", [a: "A", b: "B"])
        when:
        def encoded = encoder.encode(wfResults)
        def object = encoder.decode(encoded)

        then:
        object instanceof WorkflowOutputs
        def result = object as WorkflowOutputs
        result.createdAt == time
        result.workflowRun == "cid://1234"
        result.outputs == [a: "A", b: "B"]
    }

    def 'should encode and decode TaskRun'() {
        given:
        def encoder = new CidEncoder()
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
        def encoder = new CidEncoder()
        and:
        def time = Instant.now()
        def parameter = new Parameter("a","b", "c")
        def wfResults = new TaskOutputs("cid://1234", "cid://5678", time, [parameter], null)
        when:
        def encoded = encoder.encode(wfResults)
        def object = encoder.decode(encoded)

        then:
        object instanceof TaskOutputs
        def result = object as TaskOutputs
        result.createdAt == time
        result.taskRun == "cid://1234"
        result.workflowRun == "cid://5678"
        result.outputs.size() == 1
        result.outputs[0] == parameter
    }

    def 'object with null date attributes' () {
        given:
        def encoder = new CidEncoder()
        and:
        def wfResults = new WorkflowOutputs(null, "cid://1234")
        when:
        def encoded = encoder.encode(wfResults)
        def object = encoder.decode(encoded)
        then:
        encoded == '{"type":"WorkflowOutputs","createdAt":null,"workflowRun":"cid://1234","outputs":null,"annotations":null}'
        def result = object as WorkflowOutputs
        result.createdAt == null

    }
}
