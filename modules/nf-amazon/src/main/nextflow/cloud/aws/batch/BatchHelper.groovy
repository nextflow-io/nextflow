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

package nextflow.cloud.aws.batch

import com.amazonaws.services.batch.AWSBatch
import com.amazonaws.services.batch.model.DescribeComputeEnvironmentsRequest
import com.amazonaws.services.batch.model.DescribeJobQueuesRequest
import com.amazonaws.services.batch.model.DescribeJobsRequest
import com.amazonaws.services.ec2.AmazonEC2
import com.amazonaws.services.ec2.model.DescribeInstanceAttributeRequest
import com.amazonaws.services.ec2.model.InstanceAttributeName
import com.amazonaws.services.ecs.AmazonECS
import com.amazonaws.services.ecs.model.DescribeContainerInstancesRequest
import com.amazonaws.services.ecs.model.DescribeTasksRequest
import groovy.transform.CompileStatic
import groovy.transform.Memoized
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class BatchHelper {

    AWSBatch batchClient
    AmazonECS ecsClient
    AmazonEC2 ec2Client

    @Memoized(maxCacheSize = 100)
    protected List<String> getClusterArnByBatchQueue(String queueName) {
        final envNames = getComputeEnvByQueueName(queueName)
        return getClusterArnByCompEnvNames(envNames)
    }

    protected List<String> getClusterArnByCompEnvNames(List<String> envNames) {
        final req = new DescribeComputeEnvironmentsRequest().withComputeEnvironments(envNames)
        batchClient
                .describeComputeEnvironments(req)
                .getComputeEnvironments()
                *.getEcsClusterArn()
    }

    protected List<String> getComputeEnvByQueueName(String queueName) {
        final req = new DescribeJobQueuesRequest().withJobQueues(queueName)
        final resp = batchClient.describeJobQueues(req)
        final result = new ArrayList(10)
        for (def queue : resp.getJobQueues()) {
            for (def order : queue.getComputeEnvironmentOrder()) {
                result.add(order.getComputeEnvironment())
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
        final describeTaskReq = new DescribeTasksRequest()
                .withCluster(clusterArn)
                .withTasks(taskArn)
        final containers = ecsClient
                .describeTasks(describeTaskReq)
                .getTasks()
                *.getContainerInstanceArn()
        if( containers.size()==1 )
            return containers.get(0)
        if( containers.size()==0 )
            return null
        else
            throw new IllegalStateException("Found more than one container for taskArn=$taskArn")
    }

    protected String getInstanceIdByClusterAndContainerId(String clusterArn, String containerId) {
        final describeContainerReq = new DescribeContainerInstancesRequest()
                .withCluster(clusterArn)
                .withContainerInstances(containerId)
        final instanceIds = ecsClient
                .describeContainerInstances(describeContainerReq)
                .getContainerInstances()
                *.getEc2InstanceId()
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
        final instanceAttributeReq = new DescribeInstanceAttributeRequest()
                .withInstanceId(instanceId)
                .withAttribute(InstanceAttributeName.InstanceType)
        ec2Client
                .describeInstanceAttribute(instanceAttributeReq)
                .getInstanceAttribute()
                .getInstanceType()

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

    def describeJob(String jobId) {
        def req = new DescribeJobsRequest().withJobs(jobId)
        batchClient
                .describeJobs(req)
                .getJobs()
                .get(0)
                .getContainer()
                .getContainerInstanceArn()
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

