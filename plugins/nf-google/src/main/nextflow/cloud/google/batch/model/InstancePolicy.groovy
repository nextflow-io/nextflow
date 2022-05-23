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
 * InstancePolicy describes what instances should be created for the Job
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@ToString(includeNames = true, ignoreNulls = true, includePackage = false)
class InstancePolicy {

    // A list of allowed Compute Engine machine types, for example,
    // e2-standard-4. Default is empty which means allowing all.
    List<String> allowedMachineTypes;

    // A list of denied Compute Engine machine types.
    // Default is empty which means denying none.
    // A machine type is allowed if it matches 'allowed_machine_types' AND
    // does not match 'denied_machine_types'.
    // For example,
    //   allowed_machine_types = "e2-standard"
    //   denied_machine_types = "e2-standard-2, e2-standard-4"
    // means using all E2 standard machine types except for 'e2-standard-2' and
    // 'e2-standard-4.
    //
    // [NotImplemented]
    List<String> deniedMachineTypes;

    // A list of allowed CPU platforms, for example,
    // "Intel Cascade Lake", "AMD Rome".
    // Default is empty which means allowing all.
    //
    // [NotImplemented]
    List<String> allowedCpuPlatforms

    // A list of denied CPU platforms.
    // Default is empty which means denying none.
    // A CPU platform is allowed if it matches 'allowed_cpu_platforms' AND
    // does not match 'denied_cpu_platforms'.
    // If a CPU platform belongs to both lists, it will be denied.
    //
    // [NotImplemented]
    List<String> deniedCpuPlatforms;

    // A list of allowed accelerator types (GPU models), for example,
    // "nvidia-tesla-t4". Default is empty which means allowing all.
    //
    // [NotImplemented]
    List<String> allowedAcceleratorTypes

    // A list of denied accelerator types (GPU models).
    // Default is empty which means denying none.
    // A accelerator type is allowed if it matches 'allowed_accelerator_types'
    // AND does not match 'denied__accelerator_types'.
    //
    // [NotImplemented]
    List<String> deniedAcceleratorTypes

    // The number of accelerators per VM instance.
    //
    // [NotImplemented]
    Integer accelerator_count
    
}
