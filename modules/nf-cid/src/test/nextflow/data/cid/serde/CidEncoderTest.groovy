package nextflow.data.cid.serde

import nextflow.data.cid.model.Checksum
import nextflow.data.cid.model.DataPath
import nextflow.data.cid.model.Output
import nextflow.data.cid.model.Parameter
import nextflow.data.cid.model.TaskOutput
import nextflow.data.cid.model.TaskRun
import nextflow.data.cid.model.Workflow
import nextflow.data.cid.model.WorkflowResults
import nextflow.data.cid.model.WorkflowRun
import spock.lang.Specification

import java.time.Instant

class CidEncoderTest extends Specification{

    def 'should encode and decode Outputs'(){
        given:
            def encoder = new CidEncoder()
        and:
            def output = new TaskOutput("/path/to/file", new Checksum("hash_value", "hash_algorithm", "standard"), "cid://source", 1234)

        when:
            def encoded = encoder.encode(output)
            def object = encoder.decode(encoded)

        then:
            object instanceof Output
            def result = object as Output
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
        def workflow = new Workflow(mainScript, [], "https://nextflow.io/nf-test/", "123456")
        def wfRun = new WorkflowRun(workflow, uniqueId.toString(), "test_run", [new Parameter("String", "param1", "value1"), new Parameter("String", "param2", "value2")])

        when:
        def encoded = encoder.encode(wfRun)
        def object = encoder.decode(encoded)

        then:
        object instanceof WorkflowRun
        def result = object as WorkflowRun
        result.workflow instanceof Workflow
        result.workflow.mainScriptFile instanceof DataPath
        result.workflow.mainScriptFile.path == "file://path/to/main.nf"
        result.workflow.mainScriptFile.checksum instanceof Checksum
        result.workflow.mainScriptFile.checksum.value == "78910"
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
        def time = Instant.now().toString()
        def wfResults = new WorkflowResults(time, "cid://1234", [a: "A", b: "B"])
        when:
        def encoded = encoder.encode(wfResults)
        def object = encoder.decode(encoded)

        then:
        object instanceof WorkflowResults
        def result = object as WorkflowResults
        result.creationTime == time
        result.runId == "cid://1234"
        result.outputs == [a: "A", b: "B"]
    }

    def 'should encode and decode TaskRun'() {
        given:
        def encoder = new CidEncoder()
        and:
        def uniqueId = UUID.randomUUID()
        def taskRun = new TaskRun(
            uniqueId.toString(),"name", new Checksum("78910", "nextflow", "standard"),
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
        def time = Instant.now().toString()
        def wfResults = new WorkflowResults(time, "cid://1234", [a: "A", b: "B"])
        when:
        def encoded = encoder.encode(wfResults)
        def object = encoder.decode(encoded)

        then:
        object instanceof WorkflowResults
        def result = object as WorkflowResults
        result.creationTime == time
        result.runId == "cid://1234"
        result.outputs == [a: "A", b: "B"]
    }

    def 'object with null date attributes' () {
        given:
        def encoder = new CidEncoder()
        and:
        def wfResults = new WorkflowResults(null, "cid://1234")
        when:
        def encoded = encoder.encode(wfResults)
        def object = encoder.decode(encoded)
        then:
        encoded == '{"type":"WorkflowResults","runId":"cid://1234"}'
        def result = object as WorkflowResults
        result.creationTime == null

    }
}
