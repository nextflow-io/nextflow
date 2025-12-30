/*
 * Copyright 2020-2024, Seqera Labs
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
 *
 */

package nextflow.cloud.aws.ecs

import java.nio.file.Path
import java.util.concurrent.TimeUnit
import java.util.concurrent.TimeoutException

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.cloud.aws.AwsClientFactory
import nextflow.cloud.aws.config.AwsConfig
import nextflow.cloud.aws.nio.S3Path
import nextflow.exception.AbortOperationException
import nextflow.executor.Executor
import nextflow.fusion.FusionHelper
import nextflow.processor.ParallelPollingMonitor
import nextflow.processor.TaskHandler
import nextflow.processor.TaskMonitor
import nextflow.processor.TaskRun
import nextflow.util.Duration
import nextflow.util.RateUnit
import nextflow.util.ServiceName
import nextflow.util.ThreadPoolHelper
import nextflow.util.ThrottlingExecutor
import org.pf4j.ExtensionPoint
import software.amazon.awssdk.services.ec2.Ec2Client
import software.amazon.awssdk.services.ec2.model.DescribeSecurityGroupsRequest
import software.amazon.awssdk.services.ec2.model.DescribeSubnetsRequest
import software.amazon.awssdk.services.ec2.model.DescribeVpcsRequest
import software.amazon.awssdk.services.ec2.model.Filter
import software.amazon.awssdk.services.ecs.EcsClient
import software.amazon.awssdk.services.ecs.model.EcsException

/**
 * AWS ECS executor for running Nextflow tasks on Amazon Elastic Container Service.
 *
 * This executor uses ECS with awsvpc networking mode and requires Fusion filesystem
 * for S3-based work directory access. Tasks are executed as ECS tasks with
 * automatically registered task definitions.
 *
 * <h3>Prerequisites</h3>
 * <ul>
 *   <li>An ECS cluster (can be empty - capacity is managed automatically)</li>
 *   <li>Fusion filesystem enabled ({@code fusion.enabled = true})</li>
 *   <li>S3 bucket for work directory</li>
 *   <li>IAM execution role with permissions to pull images and write logs</li>
 *   <li>IAM task role with S3 access for the work directory</li>
 * </ul>
 *
 * <h3>Minimal IAM Execution Role Policy</h3>
 * <pre>
 * {
 *   "Version": "2012-10-17",
 *   "Statement": [
 *     {
 *       "Effect": "Allow",
 *       "Action": [
 *         "ecr:GetAuthorizationToken",
 *         "ecr:BatchCheckLayerAvailability",
 *         "ecr:GetDownloadUrlForLayer",
 *         "ecr:BatchGetImage",
 *         "logs:CreateLogStream",
 *         "logs:PutLogEvents"
 *       ],
 *       "Resource": "*"
 *     }
 *   ]
 * }
 * </pre>
 *
 * <h3>Minimal IAM Task Role Policy</h3>
 * <pre>
 * {
 *   "Version": "2012-10-17",
 *   "Statement": [
 *     {
 *       "Effect": "Allow",
 *       "Action": ["s3:GetObject", "s3:PutObject", "s3:DeleteObject", "s3:ListBucket"],
 *       "Resource": ["arn:aws:s3:::your-bucket", "arn:aws:s3:::your-bucket/*"]
 *     }
 *   ]
 * }
 * </pre>
 *
 * <h3>Minimal Configuration</h3>
 * <pre>
 * process.executor = 'awsecs'
 * fusion.enabled = true
 * wave.enabled = true
 * aws {
 *   region = 'us-east-1'
 *   ecs {
 *     cluster = 'my-cluster'
 *     executionRole = 'arn:aws:iam::ACCOUNT:role/ecsTaskExecutionRole'
 *     taskRole = 'arn:aws:iam::ACCOUNT:role/ecsTaskRole'
 *   }
 * }
 * </pre>
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@ServiceName('awsecs')
@CompileStatic
class AwsEcsExecutor extends Executor implements ExtensionPoint {

    private AwsEcsOptions awsOptions

    /**
     * ECS client instance
     */
    @PackageScope
    private EcsClient ecsClient

    /**
     * EC2 client for VPC auto-discovery
     */
    @PackageScope
    private Ec2Client ec2Client

    /**
     * AWS client factory
     */
    private AwsClientFactory clientFactory

    /**
     * Executor service to throttle service requests
     */
    private ThrottlingExecutor submitter

    /**
     * Executor service to throttle cancel requests
     */
    private ThrottlingExecutor reaper

    /**
     * Auto-discovered or configured subnets
     */
    private List<String> resolvedSubnets

    /**
     * Auto-discovered or configured security groups
     */
    private List<String> resolvedSecurityGroups

    AwsEcsOptions getAwsOptions() { awsOptions }

    @PackageScope
    EcsClient getEcsClient() { ecsClient }

    @PackageScope
    Ec2Client getEc2Client() { ec2Client }

    List<String> getResolvedSubnets() { resolvedSubnets }

    List<String> getResolvedSecurityGroups() { resolvedSecurityGroups }

    /**
     * @return {@code true} to signal containers are managed directly by AWS ECS
     */
    @Override
    final boolean isContainerNative() {
        return true
    }

    @Override
    String containerConfigEngine() {
        return 'docker'
    }

    /**
     * @return {@code true} whenever the secrets handling is managed by the executing platform itself
     */
    @Override
    final boolean isSecretNative() {
        return true
    }

    @Override
    Path getWorkDir() {
        session.bucketDir ?: session.workDir
    }

    protected void validateWorkDir() {
        /*
         * make sure the work dir is a S3 bucket
         */
        if (workDir !instanceof S3Path) {
            session.abort()
            throw new AbortOperationException("When using `$name` executor an S3 bucket must be provided as working directory using either the `-bucket-dir` or `-work-dir` command line option")
        }
    }

    protected void validatePathDir() {
        final path = session.config.navigate('env.PATH')
        if (path) {
            log.warn "Environment PATH defined in config file is ignored by AWS ECS executor"
        }
    }

    protected void validateFusion() {
        if (!isFusionEnabled()) {
            session.abort()
            throw new AbortOperationException("AWS ECS executor requires Fusion filesystem to be enabled. Please add `fusion.enabled = true` to your Nextflow configuration.")
        }
    }

    /**
     * Create AWS ECS and EC2 clients
     */
    protected void createAwsClients() {
        clientFactory = new AwsClientFactory(new AwsConfig(session.config.aws as Map))
        ecsClient = clientFactory.getEcsClient()
        ec2Client = clientFactory.getEc2Client()
        awsOptions = new AwsEcsOptions(this)
        log.debug "[AWS ECS] Executor options=$awsOptions"
    }

    /**
     * Auto-discover VPC subnets and security groups from the default VPC.
     * If no default VPC exists, fail with a clear error message.
     */
    protected void discoverVpcConfiguration() {
        // Use configured values if provided
        if (awsOptions.subnets) {
            resolvedSubnets = awsOptions.subnets
            log.debug "[AWS ECS] Using configured subnets: $resolvedSubnets"
        }
        else {
            resolvedSubnets = discoverDefaultVpcSubnets()
            log.debug "[AWS ECS] Auto-discovered subnets from default VPC: $resolvedSubnets"
        }

        if (awsOptions.securityGroups) {
            resolvedSecurityGroups = awsOptions.securityGroups
            log.debug "[AWS ECS] Using configured security groups: $resolvedSecurityGroups"
        }
        else {
            resolvedSecurityGroups = discoverDefaultVpcSecurityGroups()
            log.debug "[AWS ECS] Auto-discovered security groups from default VPC: $resolvedSecurityGroups"
        }
    }

    /**
     * Discover subnets from the default VPC.
     * @return List of subnet IDs from the default VPC
     * @throws AbortOperationException if no default VPC exists
     */
    protected List<String> discoverDefaultVpcSubnets() {
        // First, find the default VPC
        final vpcRequest = DescribeVpcsRequest.builder()
            .filters(Filter.builder().name('isDefault').values('true').build())
            .build() as DescribeVpcsRequest

        final vpcsResponse = ec2Client.describeVpcs(vpcRequest)

        if (vpcsResponse.vpcs().isEmpty()) {
            throw new AbortOperationException("""
                No default VPC found in this AWS region.

                The AWS ECS executor requires VPC configuration for task networking.
                Please configure VPC settings explicitly in your Nextflow config:

                    aws.ecs.subnets = ['subnet-xxx', 'subnet-yyy']
                    aws.ecs.securityGroups = ['sg-xxx']

                Alternatively, create a default VPC in your AWS account.
            """.stripIndent().trim())
        }

        final defaultVpcId = vpcsResponse.vpcs().first().vpcId()
        log.debug "[AWS ECS] Found default VPC: $defaultVpcId"

        // Get all subnets in the default VPC
        final subnetsRequest = DescribeSubnetsRequest.builder()
            .filters(Filter.builder().name('vpc-id').values(defaultVpcId).build())
            .build() as DescribeSubnetsRequest

        final subnetsResponse = ec2Client.describeSubnets(subnetsRequest)

        if (subnetsResponse.subnets().isEmpty()) {
            throw new AbortOperationException("""
                No subnets found in the default VPC ($defaultVpcId).

                Please configure VPC subnets explicitly in your Nextflow config:

                    aws.ecs.subnets = ['subnet-xxx', 'subnet-yyy']
            """.stripIndent().trim())
        }

        return subnetsResponse.subnets().collect { it.subnetId() }
    }

    /**
     * Discover the default security group from the default VPC.
     * @return List containing the default security group ID
     * @throws AbortOperationException if no default VPC or security group exists
     */
    protected List<String> discoverDefaultVpcSecurityGroups() {
        // First, find the default VPC
        final vpcRequest = DescribeVpcsRequest.builder()
            .filters(Filter.builder().name('isDefault').values('true').build())
            .build() as DescribeVpcsRequest

        final vpcsResponse = ec2Client.describeVpcs(vpcRequest)

        if (vpcsResponse.vpcs().isEmpty()) {
            throw new AbortOperationException("""
                No default VPC found in this AWS region.

                Please configure security groups explicitly in your Nextflow config:

                    aws.ecs.securityGroups = ['sg-xxx']
            """.stripIndent().trim())
        }

        final defaultVpcId = vpcsResponse.vpcs().first().vpcId()

        // Get the default security group for the VPC
        final sgRequest = DescribeSecurityGroupsRequest.builder()
            .filters(
                Filter.builder().name('vpc-id').values(defaultVpcId).build(),
                Filter.builder().name('group-name').values('default').build()
            )
            .build() as DescribeSecurityGroupsRequest

        final sgResponse = ec2Client.describeSecurityGroups(sgRequest)

        if (sgResponse.securityGroups().isEmpty()) {
            throw new AbortOperationException("""
                No default security group found in VPC ($defaultVpcId).

                Please configure security groups explicitly in your Nextflow config:

                    aws.ecs.securityGroups = ['sg-xxx']
            """.stripIndent().trim())
        }

        return [sgResponse.securityGroups().first().groupId()]
    }

    /**
     * Initialise the AWS ECS executor.
     */
    @Override
    protected void register() {
        super.register()
        // validate Fusion is enabled (required for ECS executor)
        validateFusion()
        // validate work directory is S3
        validateWorkDir()
        validatePathDir()
        // initialize AWS clients
        createAwsClients()
        // auto-discover VPC configuration
        discoverVpcConfiguration()
    }

    /**
     * @return The monitor instance that handles AWS ECS tasks
     */
    @Override
    protected TaskMonitor createTaskMonitor() {
        // create the throttling executor
        submitter = createExecutorService('AWSECS-executor')
        reaper = createExecutorService('AWSECS-reaper')

        final pollInterval = config.getPollInterval(name, Duration.of('10 sec'))
        final dumpInterval = config.getMonitorDumpInterval(name)
        final capacity = config.getQueueSize(name, 1000)

        final def params = [
            name: name,
            session: session,
            config: config,
            pollInterval: pollInterval,
            dumpInterval: dumpInterval,
            capacity: capacity
        ]

        log.debug "Creating parallel monitor for executor '$name' > pollInterval=$pollInterval; dumpInterval=$dumpInterval"
        new ParallelPollingMonitor(submitter, params)
    }

    /**
     * Create a task handler for the given task instance
     *
     * @param task The {@link TaskRun} instance to be executed
     * @return A {@link AwsEcsTaskHandler} for the given task
     */
    @Override
    TaskHandler createTaskHandler(TaskRun task) {
        assert task
        assert task.workDir
        log.trace "[AWS ECS] Launching process > ${task.name} -- work folder: ${task.workDirStr}"
        new AwsEcsTaskHandler(task, this)
    }

    private static final List<Integer> RETRYABLE_STATUS = [429, 500, 502, 503, 504]

    /**
     * @return Creates a {@link ThrottlingExecutor} service to throttle
     * the API requests to the AWS ECS service.
     */
    private ThrottlingExecutor createExecutorService(String name) {
        final qs = 5_000
        final limit = config.getExecConfigProp(name, 'submitRateLimit', '50/s') as String
        final size = Runtime.runtime.availableProcessors() * 5

        final opts = new ThrottlingExecutor.Options()
            .retryOn { Throwable t -> t instanceof EcsException && (t.statusCode() in RETRYABLE_STATUS) }
            .onFailure { Throwable t -> session?.abort(t) }
            .onRateLimitChange { RateUnit rate -> logRateLimitChange(rate) }
            .withRateLimit(limit)
            .withQueueSize(qs)
            .withPoolSize(size)
            .withKeepAlive(Duration.of('1 min'))
            .withAutoThrottle(true)
            .withMaxRetries(10)
            .withPoolName(name)

        ThrottlingExecutor.create(opts)
    }

    @Override
    boolean isFusionEnabled() {
        return FusionHelper.isFusionEnabled(session)
    }

    protected void logRateLimitChange(RateUnit rate) {
        log.debug "New submission rate limit: $rate"
    }

    @PackageScope
    ThrottlingExecutor getReaper() { reaper }

    @Override
    void shutdown() {
        if (submitter) {
            final tasks = submitter.shutdownNow()
            if (tasks) log.warn "Execution interrupted -- cleaning up execution pool"
            submitter.awaitTermination(5, TimeUnit.MINUTES)
        }
        // finally shutdown reaper executor
        if (reaper) {
            reaper.shutdown()
            final waitMsg = "[AWS ECS] Waiting jobs reaper to complete (%d jobs to be terminated)"
            final exitMsg = "[AWS ECS] Exiting before jobs reaper thread pool complete -- Some jobs may not be terminated"
            awaitCompletion(reaper, Duration.of('60min'), waitMsg, exitMsg)
        }
    }

    protected void awaitCompletion(ThrottlingExecutor executor, Duration duration, String waitMsg, String exitMsg) {
        try {
            ThreadPoolHelper.await(executor, duration, waitMsg, exitMsg)
        }
        catch (TimeoutException e) {
            log.warn(e.message, e)
        }
    }
}
