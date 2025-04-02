package nextflow.data.cid

import nextflow.data.cid.model.Checksum
import nextflow.data.cid.model.DataPath
import nextflow.data.cid.model.Parameter
import nextflow.data.cid.model.Workflow
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
        def configMap = [enabled: true, store: [location: storeLocation.toString(), logLocation: storeLocation.resolve(".log").toString()]]
        config = new DataConfig(configMap)
    }

    def 'should convert to Date'(){
        expect:
        CidUtils.toDate(FILE_TIME) == DATE
        where:
        FILE_TIME                   | DATE
        null                        | null
        FileTime.fromMillis(1234)   | Instant.ofEpochMilli(1234).toString()
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
        def workflow = new Workflow(mainScript, [], "https://nextflow.io/nf-test/", "123456")
        def key = "testKey"
        def value1 = new WorkflowRun(workflow, uniqueId.toString(), "test_run", [new Parameter("String", "param1", "value1"), new Parameter("String", "param2", "value2")])

        def cidStore = new DefaultCidStore()
        cidStore.open(config)
        cidStore.save(key, value1)
        when:
        List<Object> params = CidUtils.query(cidStore, new URI('cid://testKey#params'))
        then:
        params.size() == 1
        params[0] instanceof List<Parameter>
        (params[0] as List<Parameter>).size() == 2

    }

    def "should parse children elements form Fragment string"() {
        expect:
        CidUtils.parseChildrenFormFragment(FRAGMENT) == EXPECTED

        where:
        FRAGMENT            | EXPECTED
        "field1"            | ["field1"]
        "field1.field2"     | ["field1", "field2"]
        null                | []
        ""                  | []
    }

    def "should parse a query string as Map"() {
        expect:
        CidUtils.parseQuery(QUERY_STRING) == EXPECTED

        where:
        QUERY_STRING                | EXPECTED
        "key1=value1&key2=value2"   | ["key1": "value1", "key2": "value2"]
        "key=val with space"        | ["key": "val with space"]
        ""                          | [:]
        null                        | [:]
    }

    def "should check params in an object"() {
        given:
        def obj = ["field": "value", "nested": ["subfield": "subvalue"]]

        expect:
        CidUtils.checkParams(obj, PARAMS) == EXPECTED

        where:
        PARAMS                                  | EXPECTED
        ["field": "value"]                      | true
        ["field": "wrong"]                      | false
        ["nested.subfield": "subvalue"]         | true
        ["nested.subfield": "wrong"]            | false
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
        def workflow = new Workflow(mainScript, [], "https://nextflow.io/nf-test/", "123456")
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

}
