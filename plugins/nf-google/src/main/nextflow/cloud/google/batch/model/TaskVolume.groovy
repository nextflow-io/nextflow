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
 * Model a Google Batch task volume mount
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@ToString(includeNames = true, ignoreNulls = true, includePackage = false)
class TaskVolume {
    Map gcs
    String mountPath
    List<String> mountOptions

    TaskVolume withMountOptions(List<String> opts) {
        mountOptions = opts
        return this
    }

    TaskVolume withMountPath(String path) {
        this.mountPath = path
        return this
    }
}
