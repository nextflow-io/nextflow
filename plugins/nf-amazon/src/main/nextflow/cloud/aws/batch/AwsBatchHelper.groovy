/*
 * Copyright 2013-2026, Seqera Labs
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
import software.amazon.awssdk.services.ec2.model.DescribeInstancesRequest
import software.amazon.awssdk.services.ec2.model.Instance
import software.amazon.awssdk.services.ec2.model.InstanceLifecycleType
import software.amazon.awssdk.services.ecs.EcsClient
import software.amazon.awssdk.services.ecs.model.DescribeContainerInstancesRequest
import software.amazon.awssdk.services.ecs.model.DescribeTasksRequest
import software.amazon.awssdk.services.ecs.model.InvalidParameterException
import software.amazon.awssdk.services.cloudwatchlogs.CloudWatchLogsClient
import software.amazon.awssdk.services.cloudwatchlogs.model.GetLogEventsRequest
import software.amazon.awssdk.services.cloudwatchlogs.model.OutputLogEvent
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
    private BatchClient batchClient

    AwsBatchHelper(BatchClient batchClient, AwsClientFactory factory) {
        this.batchClient = batchClient
        this.factory = factory
    }

    @Memoized
    private EcsClient getEcsClient() {
        return factory.getEcsClient()
    }

    @Memoized
    private Ec2Client getEc2Client() {
        return factory.getEc2Client()
    }

    @Memoized
    private CloudWatchLogsClient getLogsClient() {
        return factory.getLogsClient()
    }

    @Memoized(maxCacheSize = 100)
    private List<String> getClusterArnByBatchQueue(String queueName) {
        final envNames = getComputeEnvByQueueName(queueName)
        return getClusterArnByCompEnvNames(envNames)
    }

    private List<String> getClusterArnByCompEnvNames(List<String> envNames) {
        final req = DescribeComputeEnvironmentsRequest.builder()
            .computeEnvironments(envNames)
            .build() as DescribeComputeEnvironmentsRequest
        batchClient
                .describeComputeEnvironments(req)
                .computeEnvironments()
                *.ecsClusterArn()
    }

    private List<String> getComputeEnvByQueueName(String queueName) {
        final req = DescribeJobQueuesRequest.builder()
            .jobQueues(queueName)
            .build() as DescribeJobQueuesRequest

        final resp = batchClient.describeJobQueues(req)

        final result = new ArrayList<String>(10)
        for (final queue : resp.jobQueues()) {
            for (final order : queue.computeEnvironmentOrder()) {
                result.add(order.computeEnvironment())
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
        final describeTaskReq = DescribeTasksRequest.builder()
            .cluster(clusterArn)
            .tasks(taskArn)
            .build() as DescribeTasksRequest
        try {
            final describeTasksResult = ecsClient.describeTasks(describeTaskReq)
            final containers =
                    describeTasksResult.tasks()
                    *.containerInstanceArn()
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
        final describeContainerReq = DescribeContainerInstancesRequest.builder()
                .cluster(clusterArn)
                .containerInstances(containerId)
                .build() as DescribeContainerInstancesRequest
        final instanceIds = ecsClient
                .describeContainerInstances(describeContainerReq)
                .containerInstances()
                *.ec2InstanceId()
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
        final req = DescribeInstancesRequest.builder()
            .instanceIds(instanceId)
            .build() as DescribeInstancesRequest
        final res = ec2Client.describeInstances(req).reservations() [0]
        final Instance instance = res ? res.instances() [0] : null
        if( !instance ) {
            log.debug "Unable to find cloud machine info for instanceId=$instanceId"
            return null
        }

        new CloudMachineInfo(
                instance.instanceType().toString(),
                instance.placement().availabilityZone(),
                getPrice(instance))
    }

    private PriceModel getPrice(Instance instance) {
        instance.instanceLifecycle() == InstanceLifecycleType.SPOT ? PriceModel.spot : PriceModel.standard
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
        final request = DescribeJobsRequest.builder()
            .jobs(jobId)
            .build() as DescribeJobsRequest
        final response = batchClient.describeJobs(request)
        if( response.jobs() ) {
            final detail = response.jobs()[0]
            return detail.container().logStreamName()
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

        final logRequest = GetLogEventsRequest.builder()
                .logGroupName(groupName ?: "/aws/batch/job")
                .logStreamName(streamId)
                .build() as GetLogEventsRequest

        final result = new StringBuilder()
        final resp = logsClient.getLogEvents(logRequest)
        for( OutputLogEvent it : resp.events() ) {
            result.append(it.message()).append('\n')
        }

        return result.toString()
    }

}

