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

    static class Accelerator {
        // The accelerator type. For example, "nvidia-tesla-t4".
        // See `gcloud compute accelerator-types list`.
        String type

        // The number of accelerators of this type.
        Integer count
    }

    static class Disk {
        // A data source from which a PD will be created.
        // Name of a public or custom image used as the data source.
        String image

        // Name of a snapshot used as the data source.
        String snapshot

        // Disk type as shown in `gcloud compute disk-types list`
        // For example, "pd-ssd", "pd-standard", "pd-balanced".
        String type

        // Disk size in GB.
        // This field is ignored if `data_source` is `disk` or `image`.
        Integer sizeGb
    }

    static class AttachedDisk {
        Disk newDisk
        // Name of an existing PD.
        String existingDisk
    }

    // The Compute Engine machine type.
    String machineType

    // The minimum CPU platform.
    // See
    // `https://cloud.google.com/compute/docs/instances/specify-min-cpu-platform`.
    String minCpuPlatform

    // The provisioning model.
    ProvisioningModel provisioningModel

    // The accelerators attached to each VM instance.
    List<Accelerator> accelerators

    // Non-boot disks to be attached for each VM created by this InstancePolicy.
    // New disks will be deleted when the attached VM is deleted.
    List<AttachedDisk> disks

    InstancePolicy withProvisioningModel(ProvisioningModel model) {
        this.provisioningModel = model
        return this
    }

    InstancePolicy withMachineType(String type) {
        this.machineType = type
        return this
    }
}
