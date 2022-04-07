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
class TaskContainer {

    /**
     * The URI to pull the container image from.
     */
    String imageUri

    /**
     * Overrides the `CMD` specified in the container. If there is an ENTRYPOINT
     * (either in the container image or with the entrypoint field below) then
     * commands are appended as arguments to the ENTRYPOINT.
     */
    List<String> commands
    /*
     * Overrides the `ENTRYPOINT` specified in the container.
     */
    String entrypoint
    /**
     * Volumes to mount (bind mount) from the host machine files or directories // into the container, formatted to match docker run's −−volume option,
     * e.g. /foo:/bar, or /foo:/bar/:ro
     */
    List<String> volumes

    /**
     * Arbitrary additional options to include in the "docker run" command when // running this container, e.g. "−−network host".
     */
    String options;

    TaskContainer withImageUri(String it) {
        this.imageUri = it
        return this
    }

    TaskContainer withEntrypoint(String it) {
        this.entrypoint = it
        return this
    }

    TaskContainer withCommands(List<String> cmd) {
        this.commands = cmd
        return this
    }

    TaskContainer withVolumes(List<String> it) {
        this.volumes = it
        return this
    }

    TaskContainer withOptions( String containerOpts ) {
        this.options = containerOpts
        return this
    }


}
