/*
 * Copyright 2022, Google Inc.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package nextflow.cloud.google.batch.model

import groovy.transform.CompileStatic
import groovy.transform.ToString

/**
 * Model Google Batch Container configuration
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@ToString(includeNames = true, ignoreNulls = true, includePackage = false)
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
