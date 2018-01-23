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

    K8sResponseException(String message, K8sResponseJson response) {
        super(message)
        this.response = response
    }

    K8sResponseException(String message, InputStream response) {
        this(message, new K8sResponseJson(fetch(response)))
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

    String getMessage() {
        response ? "${super.message}\n${response.toString().indent()}" : super.getMessage()
    }

}
