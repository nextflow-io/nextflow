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

package nextflow.cloud.aws.batch

import software.amazon.awssdk.services.batch.BatchClient
import software.amazon.awssdk.services.batch.model.DescribeComputeEnvironmentsRequest
import software.amazon.awssdk.services.batch.model.DescribeJobQueuesRequest
import software.amazon.awssdk.services.batch.model.DescribeJobsRequest
import software.amazon.awssdk.services.ec2.Ec2Client
import software.amazon.awssdk.services.ec2.model.DescribeInstanceAttributeRequest
import software.amazon.awssdk.services.ec2.model.InstanceAttributeName
import software.amazon.awssdk.services.ecs.EcsClient
import software.amazon.awssdk.services.ecs.model.DescribeContainerInstancesRequest
import software.amazon.awssdk.services.ecs.model.DescribeTasksRequest
import groovy.transform.CompileStatic
import groovy.transform.Memoized
/**
 * Implement helper method for AWS Batch
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class BatchHelper {

    BatchClient batchClient
    EcsClient ecsClient
    Ec2Client ec2Client

    @Memoized(maxCacheSize = 100)
    protected List<String> getClusterArnByBatchQueue(String queueName) {
        final envNames = getComputeEnvByQueueName(queueName)
        return getClusterArnByCompEnvNames(envNames)
    }

    protected List<String> getClusterArnByCompEnvNames(List<String> envNames) {
        final req = DescribeComputeEnvironmentsRequest.builder().computeEnvironments(envNames).build() as DescribeComputeEnvironmentsRequest
        batchClient
                .describeComputeEnvironments(req)
                .computeEnvironments()
                *.ecsClusterArn()
    }

    protected List<String> getComputeEnvByQueueName(String queueName) {
        final req = DescribeJobQueuesRequest.builder().jobQueues(queueName).build() as DescribeJobQueuesRequest
        final resp = batchClient.describeJobQueues(req)
        final result = new ArrayList(10)
        for (def queue : resp.jobQueues()) {
            for (def order : queue.computeEnvironmentOrder()) {
                result.add(order.computeEnvironment())
            }
        }
        return result
    }

    protected String getInstanceTypeByClusterAndTaskArn(String clusterArn, String taskArn) {
        final containerId = getContainerIdByClusterAndTaskArn(clusterArn, taskArn)
        final instanceId = getInstanceIdByClusterAndContainerId(clusterArn, containerId)
        instanceId ? getInstanceTypeByInstanceId(instanceId) : null
    }

    protected String getInstanceTypeByClusterAndContainerArn(String clusterArn, String containerArn) {
        final instanceId = getInstanceIdByClusterAndContainerId(clusterArn, containerArn)
        instanceId ? getInstanceTypeByInstanceId(instanceId) : null
    }

    protected String getContainerIdByClusterAndTaskArn(String clusterArn, String taskArn) {
        final describeTaskReq = DescribeTasksRequest.builder()
                .cluster(clusterArn)
                .tasks(taskArn)
                .build()
        final containers = ecsClient
                .describeTasks(describeTaskReq)
                .tasks()
                *.containerInstanceArn()
        if( containers.size()==1 )
            return containers.get(0)
        if( containers.size()==0 )
            return null
        else
            throw new IllegalStateException("Found more than one container for taskArn=$taskArn")
    }

    protected String getInstanceIdByClusterAndContainerId(String clusterArn, String containerId) {
        final describeContainerReq = DescribeContainerInstancesRequest.builder()
                .cluster(clusterArn)
                .containerInstances(containerId)
                .build()
        final instanceIds = ecsClient
                .describeContainerInstances(describeContainerReq)
                .containerInstances()
                *.ec2InstanceId()
        if( !instanceIds )
            return null
        if( instanceIds.size()==1 )
            return instanceIds.get(0)
        else
            throw new IllegalStateException("Found more than one EC2 instance for containerId=$containerId")
    }

    @Memoized(maxCacheSize = 100)
    protected String getInstanceTypeByInstanceId(String instanceId) {
        assert instanceId
        final instanceAttributeReq = DescribeInstanceAttributeRequest.builder()
                .instanceId(instanceId)
                .attribute(InstanceAttributeName.INSTANCE_TYPE)
                .build()
        ec2Client
                .describeInstanceAttribute(instanceAttributeReq)
                .instanceType()

    }

    String getInstanceTypeByQueueAndTaskArn(String queue, String taskArn) {
        final clusterArnList = getClusterArnByBatchQueue(queue)
        for( String cluster : clusterArnList ) {
            final instanceType = getInstanceTypeByClusterAndTaskArn(cluster, taskArn)
            if( instanceType )
                return instanceType
        }
        return null
    }

    String describeJob(String jobId) {
        final req = DescribeJobsRequest.builder().jobs(jobId).build()
        batchClient
                .describeJobs(req)
                .jobs()
                .get(0)
                .container()
                .containerInstanceArn()
    }

    String getInstanceTypeByQueueAndContainerArn(String queue, String containerArn) {
        final clusterArnList = getClusterArnByBatchQueue(queue)
        for( String cluster : clusterArnList ) {
            final instanceType = getInstanceTypeByClusterAndContainerArn(cluster, containerArn)
            if( instanceType )
                return instanceType
        }
        return null
    }

}

