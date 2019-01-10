/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

import groovy.transform.CompileStatic
import groovy.transform.Immutable
/**
 * Model a cloud instance
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Immutable
@CompileStatic
class CloudInstance implements Serializable, Cloneable {

    /**
     * The instance unique identifier as assigned by the provider
     */
    String id

    /**
     * The instance current state
     */
    String state

    /**
     * The instance public DNS name
     */
    String publicDnsName

    /**
     * The instance private DNS name
     */
    String privateDnsName

    /**
     * The instance private DNS name
     */
    String publicIpAddress

    /**
     * The instance private IP address
     */
    String privateIpAddress

    /**
     * The instance role in the cluster, either {@link nextflow.Const#ROLE_MASTER} or {@link nextflow.Const#ROLE_WORKER}
     */
    String role

    /**
     * The name of the cluster to which this instance belong to
     */
    String clusterName

    /**
     * @return The instance public address, falling back on the private when the public is not available
     */
    String getAddress() {
        publicDnsName ?: publicIpAddress ?: privateDnsName ?: privateIpAddress
    }

    boolean hasPublicAddress() {
        publicDnsName || publicIpAddress
    }

}
