/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
/**
 * Model a kubernetes invalid response
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@Slf4j
@CompileStatic
class K8sResponseException extends Exception {

    K8sResponseJson response

    K8sResponseException(K8sResponseJson response) {
        super(msg0(response))
    }

    K8sResponseException(String message, K8sResponseJson response) {
        super(msg1(message,response))
        this.response = response
    }

    K8sResponseException(String message, InputStream response) {
        this(message, new K8sResponseJson(fetch(response)))
    }

    static private String msg1( String msg, K8sResponseJson resp ) {
        if( !msg && resp==null )
            return null

        if( msg && resp != null ) {
            def sep = resp.isRawText() ? ' -- ' : '\n'
            return "${msg}${sep}${msg0(resp)}"
        }
        else if( msg ) {
            return msg
        }
        else {
            return msg0(resp)
        }
    }

    static private String msg0( K8sResponseJson response ) {
        if( response == null )
            return null

        if( response.isRawText() )
            response.getRawText()
        else
            "\n${response.toString().indent('  ')}"
    }

    static private String fetch(InputStream stream) {
        try {
            return stream?.text
        }
        catch( Exception e ) {
            log.debug "Unable to fetch response text -- Cause: ${e.message ?: e}"
            return null
        }
    }

}
