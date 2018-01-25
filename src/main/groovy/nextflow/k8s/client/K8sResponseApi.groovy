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

/**
 * Model a Kubernetes API response
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class K8sResponseApi {

    private int code

    private InputStream stream

    private String text

    K8sResponseApi(int code, InputStream stream) {
        this.code = code
        this.stream = stream
    }

    String toString() {
        "code=$code; stream=$stream"
    }

    int getCode() { code }

    InputStream getStream() { stream }

    String getText() {
        if( text == null ) {
            text = stream?.text
        }
        return text
    }
}
