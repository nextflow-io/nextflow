/*
 * Copyright (c) 2020-2021. Seqera Labs, S.L.
 *
 * All Rights reserved
 *
 */

package nextflow.cloud.google.batch.client

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j

/**
 * Model a Google Batch client invalid response
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@Slf4j
@CompileStatic
class BatchResponseException extends Exception {

    BatchJsonResponse response

    BatchResponseException(BatchJsonResponse response) {
        super(msg0(response))
        this.response = response
    }

    BatchResponseException(String message, BatchJsonResponse response) {
        super(msg1(message,response))
        this.response = response
    }

    BatchResponseException(String message, InputStream response) {
        this(message, new BatchJsonResponse(fetch(response)))
    }

    static private String msg1( String msg, BatchJsonResponse resp ) {
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

    static private String msg0( BatchJsonResponse response ) {
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
