import com.dnanexus.DXJSON
import com.fasterxml.jackson.databind.node.ObjectNode

/**
 * Created with IntelliJ IDEA.
 * User: bmartin
 * Date: 7/19/13
 * Time: 12:59 PM
 * To change this template use File | Settings | File Templates.
 */


ObjectNode jsonRequest = DXJSON.getObjectBuilder()
        .put("function", "process")
        .put("input", DXJSON.getObjectBuilder()
        .put("taskName", "hello")
        .put("taskEnv", "export x=1")
        .put("taskScript", "blah")
        .build())
        .build();


println jsonRequest.asText()
