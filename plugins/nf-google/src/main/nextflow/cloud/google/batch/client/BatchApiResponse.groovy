/*
 * Copyright (c) 2020-2021. Seqera Labs, S.L.
 *
 * All Rights reserved
 *
 */

package nextflow.cloud.google.batch.client

import groovy.transform.CompileStatic


/**
 * Model a Batch API response
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class BatchApiResponse {

    private int code

    private InputStream stream

    private String text

    BatchApiResponse(int code, InputStream stream) {
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
