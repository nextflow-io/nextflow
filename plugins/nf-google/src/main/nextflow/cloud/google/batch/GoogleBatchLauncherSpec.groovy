/*
 * Copyright 2023, Seqera Labs.
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

    default List<String> launchCommand() {
        return ['/bin/bash','-o','pipefail','-c', runCommand() ]
    }

    default Map<String,String> getEnvironment() {
        return Map.of()
    }
}
