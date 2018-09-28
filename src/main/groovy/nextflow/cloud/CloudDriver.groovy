/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

package nextflow.cloud
import groovy.transform.stc.ClosureParams
import groovy.transform.stc.SimpleType
import nextflow.cloud.types.CloudInstanceStatus
import nextflow.cloud.types.CloudInstanceType

import static nextflow.cloud.CloudConst.*

/**
 * Defines the cloud driver interface methods
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
trait CloudDriver {

    /**
     * Validate the cloud configuration object
     *
     * @param config The cloud launch configuration
     */
    abstract void validate( LaunchConfig config )

    /**
     * Launch cloud instances given the specified configuration
     *
     * @param instanceCount The number of instance to launch
     * @param config A {@link LaunchConfig} object providing the instance configuration
     * @return The list of launched instance IDs as given the underlying provider
     */
    abstract List<String> launchInstances( int instanceCount, LaunchConfig config )

    /**
     * Wait the list of specified instance to be in `running` state
     *
     * @param instanceIds A list of instance IDs
     * @param status The status to be reached
     */
    abstract void waitInstanceStatus( Collection<String> instanceIds, CloudInstanceStatus status )

    /**
     * Tag one or more instances with the specified key=value pairs
     *
     * @param instanceIds A list of instance IDs
     * @param tags A mpa of tags to be associated to the specified instances
     */
    abstract void tagInstances( Collection<String> instanceIds, Map<String,String> tags )

    /**
     * Tags the specified instances with the default tags
     *
     * @param instanceIds A list of instance ids
     * @param config The {@link LaunchConfig} configuration object
     */
    void tagInstances( Collection<String> instanceIds, LaunchConfig config )
    {
        final tags = [:]
        tags[TAG_CLUSTER_NAME] = config.clusterName
        tags[TAG_CLUSTER_ROLE] = config.role
        tagInstances(instanceIds, tags)
    }

    /**
     * Iterate over list of spot prices for the given instance types
     *
     * @param instanceTypes
     *      A collection of instance types
     * @param callback
     *      A closure getting a single parameter of type {@link nextflow.cloud.types.CloudSpotPrice} describing
     *      the instance properties and price for the actual instance type
     *
     */
    abstract void eachSpotPrice(
            List<String> instanceTypes,
            @ClosureParams(value=SimpleType, options = ['nextflow.cloud.types.CloudSpotPrice']) Closure callback )

    /**
     * Iterate over list of spot prices for the given instance types
     *
     * @param callback
     *      A closure getting a single parameter of type {@link nextflow.cloud.types.CloudSpotPrice} describing
     *      the instance properties and price for the actual instance type
     */
    void eachSpotPrice(
            @ClosureParams(value=SimpleType, options = ['nextflow.cloud.types.CloudSpotPrice']) Closure callback )
    {
        eachSpotPrice(Collections.emptyList(), callback)
    }

    /**
     * Iterate over the list of instances for the given tags
     *
     * @param tags
     *      One or more tag given as a {@link Map} object
     * @param callback
     *      A closure getting a single parameter of type {@link nextflow.cloud.types.CloudInstance}
     *      describing the properties for the current instance
     *
     * @see nextflow.cloud.types.CloudInstance
     */
    abstract void eachInstanceWithTags(
            Map tags,
            @ClosureParams(value=SimpleType, options = ['nextflow.cloud.types.CloudInstance']) Closure callback )

    /**
     * Iterate over the list of instances for the given instance IDs
     *
     * @param instanceIds
     *      One or more instance IDs
     * @param callback
     *      A closure getting a single parameter of type {@link nextflow.cloud.types.CloudInstance}
     *      describing the properties for the current instance
     *
     * @see nextflow.cloud.types.CloudInstance
     */
    abstract void eachInstanceWithIds(
            List<String> instanceIds,
            @ClosureParams(value=SimpleType, options = ['nextflow.cloud.types.CloudInstance']) Closure callback )

    /**
     * Iterate over the list of available instances
     *
     * @param callback
     *      A closure getting a single parameter of type {@link nextflow.cloud.types.CloudInstance}
     *      describing the properties for the current instance
     *
     * @see nextflow.cloud.types.CloudInstance
     */
    abstract void eachInstance(
            @ClosureParams(value=SimpleType, options = ['nextflow.cloud.types.CloudInstance']) Closure callback )

    /**
     * Retrieve the private IPs of all nodes in the cluster with the specified name
     *
     * @param clusterName The name of the cluster for which retrieve the list of IPs
     * @return The list of node IPs in the specified cluster
     */
    abstract List<String> listPrivateIPs( String clusterName )

    /**
     * Terminate one or more cloud instances
     *
     * @param instanceIds A list of instance IDs
     */
    abstract void terminateInstances( Collection<String> instanceIds )

    /**
     * Retrieve the instance id of the local node. This method is meant to be invoked
     * in a cloud VM
     *
     * @return The instance ID
     */
    abstract String getLocalInstanceId()

    /**
     * Retrieve a termination notice for spot/preemptive instance if available
     * This method is meant to be invoked by a cloud VM which can be retired by the
     * cloud provider
     *
     * @return
     */
    abstract String getLocalTerminationNotice()

    /**
     * Retrieve a {@link CloudInstanceType} object for the given instance type
     *
     * @param instanceType The instance type identifier e.g. {@code m4.xlarge}
     * @return The {@link CloudInstanceType} object for the given instance type
     * @throws IllegalArgumentException If the instance type is unknown
     */
    abstract CloudInstanceType describeInstanceType( String instanceType )

}
