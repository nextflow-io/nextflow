/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

package nextflow.k8s.client

import groovy.json.JsonOutput
import groovy.json.JsonSlurper
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j

/**
 * Model the response of a kubernetes api request
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class K8sResponseJson implements Map {

    @Delegate
    private Map response

    private String rawText

    K8sResponseJson(Map response) {
        this.response = response
    }

    K8sResponseJson(String response) {
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
