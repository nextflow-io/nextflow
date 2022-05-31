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
 * A Job's resource allocation policy describes when, where, and how compute
 * resources should be allocated for the Job.
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@ToString(includeNames = true, ignoreNulls = true, includePackage = false)
class AllocationPolicy {

    // Location where compute resources should be allocated for the Job.
    LocationPolicy location

    // Create only instances allowed by this policy.
    InstancePolicy instance

    // Create instances from the first instance template - MVP
    // Use 'gcloud compute instance-templates list` to see available templates
    // in the project If specified, it overrides the 'instance' field.
    Set<String> instanceTemplates

    // Create only instances in the listed provisiong models.
    // Default to allow all.
    //
    // Currently only the first model of the provisioning_models list will be
    // considered; specifying additional models (e.g., 2nd, 3rd, etc.) is a no-op.
    Set<ProvisioningModel> provisioningModels

    // Email of the service account that VMs will run as.
    String serviceAccount

    /**
     // Labels applied to all VM instances and other resources
     // created by AllocationPolicy.
     // Labels could be user provided or system generated.
     // You can assign up to 64 labels. [Google Compute Engine label
     // restrictions](https://cloud.google.com/compute/docs/labeling-resources#restrictions)
     // apply.
     // Label names that start with "goog-" or "google-" are reserved.
     */
    Map<String, String> labels

    // The network policy.
    NetworkPolicy network

    AllocationPolicy withProvisioningModel(ProvisioningModel provisioning) {
        this.provisioningModels = Collections.singleton(provisioning)
        return this
    }

    AllocationPolicy withMachineTypes(String... machineTypes) {
        this.instance = new InstancePolicy(allowedMachineTypes: machineTypes.toList())
        return this
    }

    AllocationPolicy withInstancePolicy(InstancePolicy policy) {
        this.instance = policy
        return this
    }

    AllocationPolicy withLocationPolicy(LocationPolicy policy) {
        this.location = policy
        return this
    }

    AllocationPolicy withNetworkPolicy(NetworkPolicy policy) {
        this.network = policy
        return this
    }
}
