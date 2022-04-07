/*
 * Copyright (c) 2020-2021. Seqera Labs, S.L.
 *
 * All Rights reserved
 *
 */

package nextflow.cloud.google.batch.model

import groovy.transform.CompileStatic

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class LogsPolicy {
    static public final LogsPolicy CLOUD_LOGGING = new LogsPolicy(destination: 'CLOUD_LOGGING')

    String destination
}
