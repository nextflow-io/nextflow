/*
 * Copyright 2013-2024, Seqera Labs
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

package nextflow.cloud.types

import groovy.transform.Immutable
import nextflow.util.MemoryUnit

/**
 * Models a cloud instance type
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Immutable
class CloudInstanceType implements Serializable, Cloneable {

    /**
     * Instance type ID as assigned by the cloud provider
     */
    String id

    /**
     * Number of CPUs
     */
    int cpus

    /**
     * Amount of memory (RAM)
     */
    MemoryUnit memory

    /**
     * Amount of local storage
     */
    MemoryUnit disk

    /**
     * Number of disks
     */
    int numOfDisks

    String toString() {
        "id=$id; cpus=$cpus; mem=$memory; disk=$disk x $numOfDisks"
    }
}
