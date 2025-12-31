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
import java.util.concurrent.ConcurrentHashMap

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.cloud.aws.ecs.model.ContainerDefinitionModel
import nextflow.cloud.aws.ecs.model.RegisterTaskDefinitionModel
import nextflow.exception.ProcessException
import nextflow.exception.ProcessUnrecoverableException
import nextflow.executor.BashWrapperBuilder
import nextflow.fusion.FusionAwareTask
import nextflow.processor.TaskHandler
import nextflow.processor.TaskRun
import nextflow.processor.TaskStatus
import nextflow.trace.TraceRecord
import nextflow.util.CacheHelper
import nextflow.util.MemoryUnit
import nextflow.util.TestOnly
import software.amazon.awssdk.services.ecs.EcsClient
import software.amazon.awssdk.services.ecs.model.AssignPublicIp
import software.amazon.awssdk.services.ecs.model.AwsVpcConfiguration
import software.amazon.awssdk.services.ecs.model.ContainerOverride
import software.amazon.awssdk.services.ecs.model.DescribeTasksRequest
import software.amazon.awssdk.services.ecs.model.KeyValuePair
import software.amazon.awssdk.services.ecs.model.NetworkConfiguration
import software.amazon.awssdk.services.ecs.model.RegisterTaskDefinitionResponse
import software.amazon.awssdk.services.ecs.model.RunTaskRequest
import software.amazon.awssdk.services.ecs.model.RunTaskResponse
import software.amazon.awssdk.services.ecs.model.StopTaskRequest
import software.amazon.awssdk.services.ecs.model.Task
import software.amazon.awssdk.services.ecs.model.TaskOverride

/**
 * Implements a task handler for AWS ECS Managed Instances jobs
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class AwsEcsTaskHandler extends TaskHandler implements FusionAwareTask {

    private final Path exitFile
    private final Path wrapperFile
    private final Path outputFile
    private final Path errorFile
    private final Path logFile
    private final Path scriptFile
    private final Path inputFile
    private final Path traceFile

    private AwsEcsExecutor executor

    /**
     * The ECS task ARN
     */
    private volatile String taskArn

    /**
     * The ECS task definition ARN
     */
    private volatile String taskDefArn

    /**
     * Counter for spot interruption retries
     */
    private int spotAttempts = 0

    /**
     * Cache for task definitions keyed by resource hash
     */
    private static final Map<String, String> taskDefinitionCache = new ConcurrentHashMap<>()

    @TestOnly
    protected AwsEcsTaskHandler() {}

    /**
     * Create a new ECS task handler
     *
     * @param task The {@link nextflow.processor.TaskRun} descriptor of the task to run
     * @param executor The {@link AwsEcsExecutor} instance
     */
    AwsEcsTaskHandler(TaskRun task, AwsEcsExecutor executor) {
        super(task)
        this.executor = executor
        this.logFile = task.workDir.resolve(TaskRun.CMD_LOG)
        this.scriptFile = task.workDir.resolve(TaskRun.CMD_SCRIPT)
        this.inputFile = task.workDir.resolve(TaskRun.CMD_INFILE)
        this.outputFile = task.workDir.resolve(TaskRun.CMD_OUTFILE)
        this.errorFile = task.workDir.resolve(TaskRun.CMD_ERRFILE)
        this.exitFile = task.workDir.resolve(TaskRun.CMD_EXIT)
        this.wrapperFile = task.workDir.resolve(TaskRun.CMD_RUN)
        this.traceFile = task.workDir.resolve(TaskRun.CMD_TRACE)
    }

    /**
     * @return The task ARN or null if not yet submitted
     */
    String getTaskArn() { taskArn }

    /**
     * @return The ECS client from the executor
     */
    protected EcsClient getEcsClient() { executor.ecsClient }

    /**
     * @return The AWS ECS options from the executor
     */
    protected AwsEcsOptions getAwsOptions() { executor.awsOptions }

    @Override
    boolean checkIfRunning() {
        if (taskArn == null || !isSubmitted())
            return false

        final ecsTask = describeTask(taskArn)
        if (!ecsTask)
            return false

        final lastStatus = ecsTask.lastStatus()
        final result = lastStatus in ['RUNNING', 'DEACTIVATING', 'STOPPING', 'STOPPED']

        if (result) {
            this.status = TaskStatus.RUNNING
            log.trace "[AWS ECS] Task $taskArn is now RUNNING (status: $lastStatus)"
        }

        return result
    }

    @Override
    boolean checkIfCompleted() {
        if (!isRunning())
            return false

        final ecsTask = describeTask(taskArn)
        if (!ecsTask)
            return false

        final done = ecsTask.lastStatus() == 'STOPPED'
        if (done) {
            // Extract exit code from container
            final container = ecsTask.containers()?.find { it.name() == 'main' }
            final exitCode = container?.exitCode()

            if (exitCode != null) {
                task.exitStatus = exitCode
            } else {
                // Fallback to reading exit file
                task.exitStatus = readExitFile()
            }

            // Set stdout/stderr files
            task.stdout = outputFile
            task.stderr = errorFile

            // Check for spot interruption
            final stoppedReason = ecsTask.stoppedReason()
            final stopCode = ecsTask.stopCode()?.toString()

            if (isSpotInterruption(stoppedReason, stopCode)) {
                if (spotAttempts < awsOptions.maxSpotAttempts) {
                    log.debug "[AWS ECS] Spot interruption detected for task $taskArn, attempt ${spotAttempts + 1}/${awsOptions.maxSpotAttempts}"
                    spotAttempts++
                    // Resubmit will be handled by Nextflow's retry mechanism
                }
                task.error = new ProcessException("Task terminated due to spot instance interruption: $stoppedReason")
            }
            else if (task.exitStatus != 0) {
                final reason = stoppedReason ?: 'Unknown error'
                final unrecoverable = reason.contains('CannotPullContainer') && reason.contains('unauthorized')
                task.error = unrecoverable
                    ? new ProcessUnrecoverableException(reason)
                    : new ProcessException(reason)
            }

            status = TaskStatus.COMPLETED
            log.debug "[AWS ECS] Task $taskArn completed with exit code ${task.exitStatus}"
            return true
        }

        return false
    }

    /**
     * Check if the stop reason indicates a spot instance interruption
     */
    protected boolean isSpotInterruption(String stoppedReason, String stopCode) {
        if (!stoppedReason)
            return false
        // ECS Managed Instances uses spot instances that can be interrupted
        return stoppedReason.contains('spot') ||
               stoppedReason.contains('Spot') ||
               stoppedReason.contains('SPOT') ||
               stoppedReason.contains('capacity') ||
               stopCode == 'SpotInterruption'
    }

    @Override
    void killTask() {
        if (taskArn) {
            log.trace "[AWS ECS] Stopping task: $taskArn"
            executor.reaper.submit({
                try {
                    final request = StopTaskRequest.builder()
                        .cluster(awsOptions.cluster)
                        .task(taskArn)
                        .reason('Task killed by Nextflow')
                        .build()
                    ecsClient.stopTask(request)
                    log.debug "[AWS ECS] Stopped task: $taskArn"
                }
                catch (Exception e) {
                    log.warn "[AWS ECS] Failed to stop task $taskArn: ${e.message}"
                }
            })
        }
    }

    @Override
    void submit() {
        // Build the wrapper script
        prepareLauncher()

        // Get or register task definition
        taskDefArn = getOrCreateTaskDefinition()
        log.trace "[AWS ECS] Using task definition: $taskDefArn"

        // Create and submit run task request
        final request = newRunTaskRequest()
        log.trace "[AWS ECS] Submitting task > $request"

        final response = ecsClient.runTask(request)

        if (response.tasks().isEmpty()) {
            final failures = response.failures()
            final reason = failures ? failures.collect { "${it.arn()}: ${it.reason()}" }.join('; ') : 'Unknown'
            throw new ProcessException("Failed to submit ECS task: $reason")
        }

        taskArn = response.tasks().first().taskArn()
        this.status = TaskStatus.SUBMITTED
        log.debug "[AWS ECS] Process `${task.lazyName()}` submitted > taskArn=$taskArn; work-dir=${task.workDirStr}"
    }

    @Override
    void prepareLauncher() {
        createTaskWrapper().build()
    }

    protected BashWrapperBuilder createTaskWrapper() {
        // ECS executor only supports Fusion mode
        return fusionLauncher()
    }

    /**
     * Get an existing task definition from cache or register a new one
     */
    protected String getOrCreateTaskDefinition() {
        final model = createTaskDefinitionModel()
        final cacheKey = model.computeCacheKey()

        // Check cache first
        String cachedArn = taskDefinitionCache.get(cacheKey)
        if (cachedArn) {
            log.trace "[AWS ECS] Using cached task definition: $cachedArn (key: $cacheKey)"
            return cachedArn
        }

        // Register new task definition
        final request = model.toRequest()
        log.trace "[AWS ECS] Registering task definition: ${request.family()}"

        final response = ecsClient.registerTaskDefinition(request)
        final arn = response.taskDefinition().taskDefinitionArn()

        // Cache the result
        taskDefinitionCache.put(cacheKey, arn)
        log.debug "[AWS ECS] Registered task definition: $arn (key: $cacheKey)"

        return arn
    }

    /**
     * Create a task definition model from the current task configuration
     */
    protected RegisterTaskDefinitionModel createTaskDefinitionModel() {
        final model = RegisterTaskDefinitionModel.create()

        // Set task definition family name
        final containerImage = task.container
        final familyName = "nf-${sanitizeFamilyName(containerImage)}"
        model.family = familyName

        // Set IAM roles
        model.executionRoleArn = awsOptions.executionRole
        model.taskRoleArn = awsOptions.taskRole

        // Map CPU (Nextflow cpus -> ECS CPU units, 1 CPU = 1024 units)
        final cpus = task.config.getCpus()
        model.cpu = (cpus * 1024).toString()

        // Map memory (Nextflow memory -> ECS memory in MiB)
        final memory = task.config.getMemory()
        model.memory = memory ? (memory.toMega() as Long).toString() : '1024'

        // Set ephemeral storage if specified
        final disk = task.config.getDisk()
        if (disk) {
            model.ephemeralStorageGiB = disk.toGiga() as Integer
        }

        // Configure container
        final container = model.containerDefinitions.first()
        container.image = containerImage
        container.withLogging(awsOptions.logsGroup, 'nf', awsOptions.region)

        // Set container resource limits (soft limits at container level)
        container.cpu = cpus * 1024
        container.memory = memory ? memory.toMega() as Integer : 1024

        // Enable Fusion filesystem support (adds SYS_ADMIN capability and /dev/fuse device)
        // This is always enabled because ECS executor requires Fusion
        container.fusionEnabled = true

        return model
    }

    /**
     * Sanitize container image name for use in task definition family name
     */
    protected String sanitizeFamilyName(String imageName) {
        if (!imageName)
            return 'unknown'
        // Replace invalid characters with hyphens, keep only alphanumeric and hyphens
        def result = imageName
            .replaceAll('[/:@]', '-')
            .replaceAll('[^a-zA-Z0-9-]', '')
            .replaceAll('-+', '-')
            .replaceAll('^-|-$', '')
        // Limit length (ECS family names max 255 chars)
        if (result.length() > 200)
            result = result.substring(0, 200)
        return result ?: 'task'
    }

    /**
     * Create the RunTask request
     */
    protected RunTaskRequest newRunTaskRequest() {
        // Build network configuration
        final assignPublicIp = awsOptions.assignPublicIp
            ? AssignPublicIp.ENABLED
            : AssignPublicIp.DISABLED

        final networkConfig = NetworkConfiguration.builder()
            .awsvpcConfiguration(
                AwsVpcConfiguration.builder()
                    .subnets(executor.resolvedSubnets)
                    .securityGroups(executor.resolvedSecurityGroups)
                    .assignPublicIp(assignPublicIp)
                    .build()
            )
            .build()

        // Build container overrides with command and environment
        final envVars = getEnvironment().collect { k, v ->
            KeyValuePair.builder().name(k).value(v).build()
        }

        final containerOverride = ContainerOverride.builder()
            .name('main')
            .command(getContainerCommand())
            .environment(envVars)
            .build()

        final taskOverride = TaskOverride.builder()
            .containerOverrides(containerOverride)
            .build()

        return RunTaskRequest.builder()
            .cluster(awsOptions.cluster)
            .taskDefinition(taskDefArn)
            .networkConfiguration(networkConfig)
            .overrides(taskOverride)
            .count(1)
            .build()
    }

    /**
     * Get the command to execute in the container
     */
    protected List<String> getContainerCommand() {
        // ECS executor only supports Fusion mode
        return fusionSubmitCli()
    }

    /**
     * Get environment variables for the container
     */
    protected Map<String, String> getEnvironment() {
        final result = new LinkedHashMap<String, String>()

        // Add Fusion environment if enabled
        if (fusionEnabled()) {
            for (Map.Entry<String, String> entry : fusionLauncher().fusionEnv()) {
                result.put(entry.key, entry.value)
            }
        }

        return result
    }

    /**
     * Describe an ECS task by ARN
     */
    protected Task describeTask(String taskArn) {
        try {
            final request = DescribeTasksRequest.builder()
                .cluster(awsOptions.cluster)
                .tasks(taskArn)
                .build()

            final response = ecsClient.describeTasks(request)

            if (response.tasks().isEmpty()) {
                log.debug "[AWS ECS] Task not found: $taskArn"
                return null
            }

            return response.tasks().first()
        }
        catch (Exception e) {
            log.warn "[AWS ECS] Failed to describe task $taskArn: ${e.message}"
            return null
        }
    }

    /**
     * Read exit code from the exit file
     */
    protected int readExitFile() {
        try {
            exitFile.text as Integer
        }
        catch (Exception e) {
            log.debug "[AWS ECS] Cannot read exit status for task: `${task.lazyName()}` | ${e.message}"
            return Integer.MAX_VALUE
        }
    }

    @Override
    TraceRecord getTraceRecord() {
        def result = super.getTraceRecord()
        result.put('native_id', taskArn)
        return result
    }
}
