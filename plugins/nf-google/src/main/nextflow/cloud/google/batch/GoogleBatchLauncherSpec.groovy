/*
 * Copyright (c) 2020-2021. Seqera Labs, S.L.
 *
 * All Rights reserved
 *
 */

package nextflow.cloud.google.batch

import com.google.cloud.batch.v1.Volume

/**
 * Defines the operation supported by Google Batch tasks launcher
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
interface GoogleBatchLauncherSpec {

    /**
     * @return
     *      A list of string representing the container mounts. Each mount uses the docker
     *      mount conventional syntax e.g. {@code /mnt/disks/foo/scratch:/mnt/disks/foo/scratch:rw}
     */
    List<String> getContainerMounts()

    /**
     * @return A list of Batch {@link Volume} to be made accessible to the container
     */
    List<Volume> getVolumes()

    /**
     * @return A string representing the command to be executed by the containerised task
     */
    String runCommand()
}
