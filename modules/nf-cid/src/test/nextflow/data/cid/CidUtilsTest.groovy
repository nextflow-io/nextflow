package nextflow.data.cid

import nextflow.data.cid.model.Checksum
import nextflow.data.cid.model.DataPath
import nextflow.data.cid.model.Parameter
import nextflow.data.cid.model.Workflow
import nextflow.data.cid.model.WorkflowOutputs
import nextflow.data.cid.model.WorkflowRun
import nextflow.data.config.DataConfig
import spock.lang.Specification
import spock.lang.TempDir

import java.nio.file.Path
import java.nio.file.attribute.FileTime
import java.time.Instant

class CidUtilsTest extends Specification{

    @TempDir
    Path tempDir

    Path storeLocation
    DataConfig config

    def setup() {
        storeLocation = tempDir.resolve("store")
        def configMap = [enabled: true, store: [location: storeLocation.toString()]]
        config = new DataConfig(configMap)
    }

    def 'should convert to Date'(){
        expect:
        CidUtils.toDate(FILE_TIME) == DATE
        where:
        FILE_TIME                   | DATE
        null                        | null
        FileTime.fromMillis(1234)   | Instant.ofEpochMilli(1234)
    }

    def 'should convert to FileTime'(){
        expect:
        CidUtils.toFileTime(DATE) == FILE_TIME
        where:
        FILE_TIME                   | DATE
        null                        | null
        FileTime.fromMillis(1234)   | Instant.ofEpochMilli(1234).toString()
    }


    def 'should query'() {
        given:
        def uniqueId = UUID.randomUUID()
        def mainScript = new DataPath("file://path/to/main.nf", new Checksum("78910", "nextflow", "standard"))
        def workflow = new Workflow([mainScript], "https://nextflow.io/nf-test/", "123456")
        def key = "testKey"
        def value1 = new WorkflowRun(workflow, uniqueId.toString(), "test_run", [new Parameter("String", "param1", "value1"), new Parameter("String", "param2", "value2")])
        def outputs1 = new WorkflowOutputs(Instant.now(), "cid://testKey", [output: "name"] )
        def cidStore = new DefaultCidStore()
        cidStore.open(config)
        cidStore.save(key, value1)
        cidStore.save("$key/outputs", outputs1)

        when:
        List<Object> params = CidUtils.query(cidStore, new URI('cid://testKey#params'))
        then:
        params.size() == 1
        params[0] instanceof List<Parameter>
        (params[0] as List<Parameter>).size() == 2

        when:
        List<Object> outputs = CidUtils.query(cidStore, new URI('cid://testKey#outputs'))
        then:
        outputs.size() == 1
        outputs[0] instanceof Map
        outputs[0]['output'] == "name"

        expect:
        CidUtils.query(cidStore, new URI('cid://testKey#no-exist')) == []
        CidUtils.query(cidStore, new URI('cid://testKey#outputs.no-exist')) == []
        CidUtils.query(cidStore, new URI('cid://no-exist#something')) == []

    }

    def "should parse children elements form Fragment string"() {
        expect:
        CidUtils.parseChildrenFormFragment(FRAGMENT) == EXPECTED as String[]

        where:
        FRAGMENT                | EXPECTED
        "workflow"              | ["workflow"]
        "workflow.repository"   | ["workflow", "repository"]
        null                    | []
        ""                      | []
    }

    def "should parse a query string as Map"() {
        expect:
        CidUtils.parseQuery(QUERY_STRING) == EXPECTED

        where:
        QUERY_STRING                | EXPECTED
        "type=value1&taskRun=value2"   | ["type": "value1", "taskRun": "value2"]
        "type=val with space"        | ["type": "val with space"]
        ""                          | [:]
        null                        | [:]
    }

    def "should check params in an object"() {
        given:
        def obj = [ "type": "value", "workflow": ["repository": "subvalue"], "outputs" : [ ["path":"/to/file"],["path":"file2"] ] ]

        expect:
        CidUtils.checkParams(obj, PARAMS) == EXPECTED

        where:
        PARAMS                                  | EXPECTED
        ["type": "value"]                       | true
        ["type": "wrong"]                       | false
        ["workflow.repository": "subvalue"]     | true
        ["workflow.repository": "wrong"]        | false
        ["outputs.path": "wrong"]               | false
        ["outputs.path": "/to/file"]            | true
        ["outputs.path": "file2"]               | true

    }

    def 'should parse query' (){
        expect:
        CidUtils.parseQuery(PARAMS) == EXPECTED
        where:
        PARAMS                              | EXPECTED
        "type=value"                        | ["type": "value"]
        "workflow.repository=subvalue"      | ["workflow.repository": "subvalue"]
        ""                                  | [:]
        null                                | [:]
    }

    def "should navigate in object params"() {
        given:
        def obj = [
            "key1": "value1",
            "nested": [
                "subkey": "subvalue"
            ]
        ]

        expect:
        CidUtils.navigate(obj, PATH) == EXPECTED

        where:
        PATH             | EXPECTED
        "key1"           | "value1"
        "nested.subkey"  | "subvalue"
        "wrongKey"       | null
    }

    def "should add objects matching parameters"() {
        given:
        def results = []

        when:
        CidUtils.treatObject(OBJECT, PARAMS, results)

        then:
        results == EXPECTED

        where:
        OBJECT                                                                  | PARAMS                            | EXPECTED
        ["field": "value"]                                                      | ["field": "value"]                | [["field": "value"]]
        ["field": "wrong"]                                                      | ["field": "value"]                | []
        [["field": "value"], ["field": "x"]]                                    | ["field": "value"]                | [["field": "value"]]
        "string"                                                                | [:]                               | ["string"]
        ["nested": ["subfield": "match"]]                                       | ["nested.subfield": "match"]      | [["nested": ["subfield": "match"]]]
        ["nested": ["subfield": "nomatch"]]                                     | ["nested.subfield": "match"]      | []
        [["nested": ["subfield": "match"]], ["nested": ["subfield": "other"]]]  | ["nested.subfield": "match"]      | [["nested": ["subfield": "match"]]]
    }

    def "Should search path"() {
        given:
        def uniqueId = UUID.randomUUID()
        def mainScript = new DataPath("file://path/to/main.nf", new Checksum("78910", "nextflow", "standard"))
        def workflow = new Workflow([mainScript], "https://nextflow.io/nf-test/", "123456")
        def key = "testKey"
        def value1 = new WorkflowRun(workflow, uniqueId.toString(), "test_run", [new Parameter("String", "param1", "value1"), new Parameter("String", "param2", "value2")])
        def cidStore = new DefaultCidStore()
        cidStore.open(config)
        cidStore.save(key, value1)
        when:
        def result = CidUtils.searchPath(cidStore, key, ["name":"param1"], ["params"] as String[])

        then:
        result == [new Parameter("String", "param1", "value1")]
    }

    def 'should navigate' (){
        def uniqueId = UUID.randomUUID()
        def mainScript = new DataPath("file://path/to/main.nf", new Checksum("78910", "nextflow", "standard"))
        def workflow = new Workflow([mainScript], "https://nextflow.io/nf-test/", "123456")
        def wfRun = new WorkflowRun(workflow, uniqueId.toString(), "test_run", [new Parameter("String", "param1", [key: "value1"]), new Parameter("String", "param2", "value2")])

        expect:
            CidUtils.navigate(wfRun, "workflow.commitId") == "123456"
            CidUtils.navigate(wfRun, "params.name") == ["param1", "param2"]
            CidUtils.navigate(wfRun, "params.value.key") == "value1"
            CidUtils.navigate(wfRun, "params.value.no-exist") == null
            CidUtils.navigate(wfRun, "params.no-exist") == null
            CidUtils.navigate(wfRun, "no-exist") == null
            CidUtils.navigate(null, "something") == null
    }

}
