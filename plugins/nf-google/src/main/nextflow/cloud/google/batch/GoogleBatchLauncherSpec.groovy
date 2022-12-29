/*
 * Copyright (c) 2020-2021. Seqera Labs, S.L.
 *
 * All Rights reserved
 *
 */

package nextflow.cloud.google.batch

import com.google.cloud.batch.v1.Volume

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
interface GoogleBatchLauncherSpec {

    List<String> getContainerMounts()

    List<Volume> getVolumes()

    String runCommand()
}
