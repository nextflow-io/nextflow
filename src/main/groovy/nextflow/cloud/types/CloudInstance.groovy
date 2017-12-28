/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
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
