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

    static class InstancePolicyOrTemplate {
        InstancePolicy policy
        // Name of an instance template used to create VMs.
        // Named the field as 'instance_template' instead of 'template' to avoid
        // c++ keyword conflict.
        String instance_template

        InstancePolicyOrTemplate(InstancePolicy policy) {
            this.policy = policy
        }
    }

    // Location where compute resources should be allocated for the Job.
    LocationPolicy location

    // Describe instances that can be created by this AllocationPolicy.
    // Only instances[0] is supported now.
    List<InstancePolicyOrTemplate> instance

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

    InstancePolicy instancePolicy() {
        if( !instance ) {
            this.instance = List.of(new InstancePolicyOrTemplate(new InstancePolicy()))
        }
        return instance[0].policy
    }

    AllocationPolicy withProvisioningModel(ProvisioningModel provisioning) {
        instancePolicy().withProvisioningModel(provisioning)
        return this
    }

    AllocationPolicy withMachineType(String machineType) {
        instancePolicy().withMachineType(machineType)
        return this
    }

    AllocationPolicy withNetworkPolicy(NetworkPolicy policy) {
        this.network = policy
        return this
    }
}
