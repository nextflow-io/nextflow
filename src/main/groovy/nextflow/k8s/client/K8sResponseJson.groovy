/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
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
