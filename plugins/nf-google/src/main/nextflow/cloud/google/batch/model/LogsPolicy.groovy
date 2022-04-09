/*
 * Copyright (c) 2020-2021. Seqera Labs, S.L.
 *
 * All Rights reserved
 *
 */

package nextflow.cloud.google.batch.model

import groovy.transform.CompileStatic
import groovy.transform.ToString

/**
 * Model Google Batch Logs policy
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@ToString(includeNames = true, ignoreNulls = true, includePackage = false)
class LogsPolicy {
    static public final LogsPolicy CLOUD_LOGGING = new LogsPolicy(destination: 'CLOUD_LOGGING')

    String destination
}
