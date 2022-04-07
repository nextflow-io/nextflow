/*
 * Copyright (c) 2020-2021. Seqera Labs, S.L.
 *
 * All Rights reserved
 *
 */

package nextflow.cloud.google.batch.client

import groovy.json.JsonOutput
import groovy.json.JsonSlurper
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j

/**
 * Model the response of a Google Batch api request
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class BatchJsonResponse implements Map {

    @Delegate
    private Map response

    private String rawText

    BatchJsonResponse(Map response) {
        this.response = response
    }

    BatchJsonResponse(String response) {
        this.response = toJson(response)
        this.rawText = response
    }

    boolean isRawText() { !response && rawText }

    String getRawText() { rawText }

    static private Map toJson(String raw) {
        try {
            return (Map)new JsonSlurper().parseText(raw)
        }
        catch( Exception e ) {
            log.trace "[K8s] cannot parse response to json -- raw: ${raw? '\n'+raw.indent('  ') :'null'}"
            return Collections.emptyMap()
        }
    }

    static private String prettyPrint(String json) {
        try {
            JsonOutput.prettyPrint(json)
        }
        catch( Exception e ) {
            return json
        }
    }

    String toString() {
        response ? prettyPrint(JsonOutput.toJson(response)) : rawText
    }

}
