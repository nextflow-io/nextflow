package nextflow.data.cid

import nextflow.data.cid.model.Checksum
import nextflow.data.cid.model.DataPath
import nextflow.data.cid.model.Parameter
import nextflow.data.cid.model.Workflow
import nextflow.data.cid.model.WorkflowOutput
import nextflow.data.cid.model.WorkflowRun
import nextflow.data.config.DataConfig
import spock.lang.Specification
import spock.lang.TempDir

import java.nio.file.Path

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
        List<Object> params = CidUtils.query(cidStore, new URI('cid://testKey/params'))
        then:
        params.size() == 1
        params[0] instanceof List<Parameter>
        (params[0] as List<Parameter>).size() == 2

    }

}
