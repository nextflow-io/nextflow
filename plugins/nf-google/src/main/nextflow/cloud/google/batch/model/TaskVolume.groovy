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
