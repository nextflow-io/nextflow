/*
 * Copyright 2021, Microsoft Corp
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
package nextflow.cloud.azure.batch

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString
import nextflow.util.MemoryUnit

/**
 * Model the size of a Azure VM
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@ToString(includeNames = true, includePackage = false)
@EqualsAndHashCode
@CompileStatic
class AzVmType {
    String name
    Integer maxDataDiskCount
    MemoryUnit memory
    Integer numberOfCores
    MemoryUnit osDiskSize
    MemoryUnit resourceDiskSize

    AzVmType() {}

    AzVmType(Map map) {
        this.name = map.name
        this.maxDataDiskCount = map.maxDataDiskCount as Integer
        this.memory = map.memoryInMB ? MemoryUnit.of( "$map.memoryInMB MB" ) : null
        this.numberOfCores = map.numberOfCores as Integer
        this.osDiskSize = map.osDiskSizeInMB ? MemoryUnit.of( "$map.osDiskSizeInMB MB" ) : null
        this.resourceDiskSize = map.resourceDiskSizeInMB ? MemoryUnit.of( "$map.resourceDiskSizeInMB MB" ) : null
    }
}
