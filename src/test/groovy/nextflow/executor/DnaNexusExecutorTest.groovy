package nextflow.executor
import com.dnanexus.DXJSON
import com.fasterxml.jackson.databind.JsonNode
import com.fasterxml.jackson.databind.node.ObjectNode
import nextflow.util.DxHelper
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DnaNexusExecutorTest extends Specification {



    def testCreateInputObject() {

        setup:
        ObjectNode processJobInputHash

        String instance = "dx_m1.large"
        JsonNode input = DxHelper.objToJson([ makeDXLink('file-xxx') ])
        JsonNode output = DxHelper.objToJson([])

        String taskInputId =  "file-t4sk1np0t"
        String scriptId = "file-B8fx3b801xqvy0gbP9k084F9"
        String name = "task1(1)"
        String env = " export CLASSPATH='/usr/share/java/dnanexus-api-0.1.0.jar:' export DNANEXUS_HOME='/usr/share/dnanexus' "

        if (taskInputId){
            processJobInputHash = DXJSON.getObjectBuilder()
                 .put("function", "process")
                 .put("input", DXJSON.getObjectBuilder()
                     .put("inputs", input)
                     .put("outputs", output)
                     .put("taskName", name)
                     .put("taskEnv", env)
                     .put("taskScript", makeDXLink(scriptId))
                     .put("taskInput", makeDXLink(taskInputId))
                     .build())
                 .put("systemRequirements", DXJSON.getObjectBuilder()
                     .put("process", DXJSON.getObjectBuilder()          // "*"
                         .put("instanceType", instance)
                         .build())
                     .build())
                 .build()
        }
        else{
            processJobInputHash = DXJSON.getObjectBuilder()
                 .put("function", "process")
                 .put("input", DXJSON.getObjectBuilder()
                     .put("inputs", input)
                     .put("outputs", output)
                     .put("taskName", name)
                     .put("taskEnv", env)
                     .put("taskScript", makeDXLink(scriptId))
                     .build())
                 .put("systemRequirements", DXJSON.getObjectBuilder()
                     .put("process", DXJSON.getObjectBuilder()      // "*"
                         .put("instanceType", instance)
                         .build())
                     .build())
                 .build()
        }

        def expectedJson = processJobInputHash


        when:
        def myJson = DnaNexusExecutor.createInputObject(input.toList(),output.toList(),name,env,scriptId,taskInputId,instance)

        then:
        myJson == expectedJson

    }

    protected ObjectNode makeDXLink(String objectId) {
        return DXJSON.getObjectBuilder().put('$dnanexus_link', objectId).build();
    }

}
