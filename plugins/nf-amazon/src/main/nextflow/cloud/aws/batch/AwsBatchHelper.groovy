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

import com.amazonaws.services.batch.AWSBatch
import com.amazonaws.services.batch.model.DescribeComputeEnvironmentsRequest
import com.amazonaws.services.batch.model.DescribeJobQueuesRequest
import com.amazonaws.services.batch.model.DescribeJobsRequest
import com.amazonaws.services.ec2.AmazonEC2
import com.amazonaws.services.ec2.model.DescribeInstancesRequest
import com.amazonaws.services.ec2.model.Instance
import com.amazonaws.services.ecs.AmazonECS
import com.amazonaws.services.ecs.model.DescribeContainerInstancesRequest
import com.amazonaws.services.ecs.model.DescribeTasksRequest
import com.amazonaws.services.ecs.model.InvalidParameterException
import com.amazonaws.services.logs.AWSLogs
import com.amazonaws.services.logs.model.GetLogEventsRequest
import com.amazonaws.services.logs.model.OutputLogEvent
import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.util.logging.Slf4j
import nextflow.cloud.aws.AwsClientFactory
import nextflow.cloud.types.CloudMachineInfo
import nextflow.cloud.types.PriceModel
/**
 * Helper class to resolve Batch related metadata
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class AwsBatchHelper {

    private AwsClientFactory factory
    private AWSBatch batchClient

    AwsBatchHelper(AWSBatch batchClient, AwsClientFactory factory) {
        this.batchClient = batchClient
        this.factory = factory
    }

    @Memoized
    private AmazonECS getEcsClient() {
        return factory.getEcsClient()
    }

    @Memoized
    private AmazonEC2 getEc2Client() {
        return factory.getEc2Client()
    }

    @Memoized
    private AWSLogs getLogsClient() {
        return factory.getLogsClient()
    }

    @Memoized(maxCacheSize = 100)
    private List<String> getClusterArnByBatchQueue(String queueName) {
        final envNames = getComputeEnvByQueueName(queueName)
        return getClusterArnByCompEnvNames(envNames)
    }

    private List<String> getClusterArnByCompEnvNames(List<String> envNames) {
        final req = new DescribeComputeEnvironmentsRequest().withComputeEnvironments(envNames)
        batchClient
                .describeComputeEnvironments(req)
                .getComputeEnvironments()
                *.getEcsClusterArn()
    }

    private List<String> getComputeEnvByQueueName(String queueName) {
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

    private CloudMachineInfo getInfoByClusterAndTaskArn(String clusterArn, String taskArn) {
        final containerId = getContainerIdByClusterAndTaskArn(clusterArn, taskArn)
        final instanceId = containerId ? getInstanceIdByClusterAndContainerId(clusterArn, containerId) : null as String
        return instanceId ? getInfoByInstanceId(instanceId) : null
    }

    private String getContainerIdByClusterAndTaskArn(String clusterArn, String taskArn) {
        final describeTaskReq = new DescribeTasksRequest()
                .withCluster(clusterArn)
                .withTasks(taskArn)
        try {
            final describeTasksResult = ecsClient.describeTasks(describeTaskReq)
            final containers =
                    describeTasksResult.getTasks()
                    *.getContainerInstanceArn()
            if( containers.size()==1 ) {
                return containers.get(0)
            }
            if( containers.size()==0 ) {
                log.debug "Unable to find container id for clusterArn=$clusterArn and taskArn=$taskArn"
                return null
            }
            else
                throw new IllegalStateException("Found more than one container for taskArn=$taskArn")
        }
        catch (InvalidParameterException e) {
            log.debug "Cannot find container id for clusterArn=$clusterArn and taskArn=$taskArn - The task is likely running on another cluster"
            return null
        }
    }

    private String getInstanceIdByClusterAndContainerId(String clusterArn, String containerId) {
        final describeContainerReq = new DescribeContainerInstancesRequest()
                .withCluster(clusterArn)
                .withContainerInstances(containerId)
        final instanceIds = ecsClient
                .describeContainerInstances(describeContainerReq)
                .getContainerInstances()
                *.getEc2InstanceId()
        if( !instanceIds ) {
            log.debug "Unable to find EC2 instance id for clusterArn=$clusterArn and containerId=$containerId"
            return null
        }
        if( instanceIds.size()==1 )
            return instanceIds.get(0)
        else
            throw new IllegalStateException("Found more than one EC2 instance for containerId=$containerId")
    }

    @Memoized(maxCacheSize = 1_000)
    private CloudMachineInfo getInfoByInstanceId(String instanceId) {
        assert instanceId
        final req = new DescribeInstancesRequest() .withInstanceIds(instanceId)
        final res = ec2Client .describeInstances(req) .getReservations() [0]
        final Instance instance = res ? res.getInstances() [0] : null
        if( !instance ) {
            log.debug "Unable to find cloud machine info for instanceId=$instanceId"
            return null
        }

        new CloudMachineInfo(
                instance.getInstanceType(),
                instance.getPlacement().getAvailabilityZone(),
                getPrice(instance))
    }

    private PriceModel getPrice(Instance instance) {
        instance.getInstanceLifecycle()=='spot' ? PriceModel.spot : PriceModel.standard
    }

    CloudMachineInfo getCloudInfoByQueueAndTaskArn(String queue, String taskArn) {
        final clusterArnList = getClusterArnByBatchQueue(queue)
        for( String cluster : clusterArnList ) {
            final result = getInfoByClusterAndTaskArn(cluster, taskArn)
            if( result )
                return result
        }

        log.debug "Unable to find cloud info for queue=$queue and taskArn=$taskArn"
        return null
    }

    protected String getLogStreamId(String jobId) {
        final request = new DescribeJobsRequest() .withJobs(jobId)
        final response = batchClient.describeJobs(request)
        if( response.jobs ) {
            final detail = response.jobs[0]
            return detail.container.logStreamName
        }
        else {
            log.debug "Unable to find info for batch job id=$jobId"
            return null
        }
    }

    /**
     * Retrieve the cloudwatch logs for the specified AWS Batch Job ID
     *
     * @param jobId
     *      The Batch Job ID for which retrieve the job
     * @return
     *      The Batch jobs as a string value or {@code null} if no logs is available. Note, if the log
     *      is made of multiple *page* this method returns only the first one
     */
    String getTaskLogStream(String jobId, String groupName) {
        final streamId = getLogStreamId(jobId)
        if( !streamId ) {
            log.debug "Unable to find CloudWatch log stream for batch job id=$jobId"
            return null
        }

        final logRequest = new GetLogEventsRequest()
                .withLogGroupName(groupName ?: "/aws/batch/job")
                .withLogStreamName(streamId)

        final result = new StringBuilder()
        final resp = logsClient .getLogEvents(logRequest)
        for( OutputLogEvent it : resp.events ) {
            result.append(it.getMessage()).append('\n')
        }
        
        return result.toString()
    }

}

