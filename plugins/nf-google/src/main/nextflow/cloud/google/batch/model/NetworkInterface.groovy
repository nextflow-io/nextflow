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
 * Model a network interface
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@ToString(includeNames = true, ignoreNulls = true, includePackage = false)
class NetworkInterface {

    /*
     * The URL of the network resource.
     */
    String network

    /**
     * The URL of the Subnetwork resource.
     */
    String subnetwork

    /**
     // Default is false (with an external IP address). Required if
     // no external public IP address is attached to the VM. If no external
     // public IP address, additional configuration is required to allow the VM
     // to access Google Services. See
     // https://cloud.google.com/vpc/docs/configure-private-google-access and
     // https://cloud.google.com/nat/docs/gce-example#create-nat for more
     // information.
     */
    Boolean noExternalIpAddress

    NetworkInterface withNetwork(String network) {
        this.network = network
        return this
    }

    NetworkInterface withSubnetwork(String subnetwork) {
        this.subnetwork = subnetwork
        return this
    }

    NetworkInterface withNoExternalIpAddress(boolean value) {
        this.noExternalIpAddress = value
        return this
    }
}
