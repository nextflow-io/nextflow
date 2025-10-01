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


import static nextflow.cloud.aws.batch.AwsContainerOptionsMapper.*

import java.nio.file.Path
import java.nio.file.Paths
import java.time.Instant

import groovy.transform.Canonical
import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.util.logging.Slf4j
import nextflow.BuildInfo
import nextflow.SysEnv
import nextflow.cloud.aws.batch.model.ContainerPropertiesModel
import nextflow.cloud.aws.batch.model.RegisterJobDefinitionModel
import nextflow.cloud.types.CloudMachineInfo
import nextflow.container.ContainerNameValidator
import nextflow.exception.ProcessException
import nextflow.exception.ProcessSubmitException
import nextflow.exception.ProcessUnrecoverableException
import nextflow.executor.BashWrapperBuilder
import nextflow.fusion.FusionAwareTask
import nextflow.processor.BatchContext
import nextflow.processor.BatchHandler
import nextflow.processor.TaskArrayRun
import nextflow.processor.TaskHandler
import nextflow.processor.TaskRun
import nextflow.processor.TaskStatus
import nextflow.trace.TraceRecord
import nextflow.util.CacheHelper
import nextflow.util.MemoryUnit
import nextflow.util.TestOnly
import software.amazon.awssdk.services.batch.BatchClient
import software.amazon.awssdk.services.batch.model.ArrayProperties
import software.amazon.awssdk.services.batch.model.AssignPublicIp
import software.amazon.awssdk.services.batch.model.AttemptContainerDetail
import software.amazon.awssdk.services.batch.model.BatchException
import software.amazon.awssdk.services.batch.model.ClientException
import software.amazon.awssdk.services.batch.model.ContainerOverrides
import software.amazon.awssdk.services.batch.model.DescribeJobDefinitionsRequest
import software.amazon.awssdk.services.batch.model.DescribeJobDefinitionsResponse
import software.amazon.awssdk.services.batch.model.DescribeJobsRequest
import software.amazon.awssdk.services.batch.model.EphemeralStorage
import software.amazon.awssdk.services.batch.model.EvaluateOnExit
import software.amazon.awssdk.services.batch.model.Host
import software.amazon.awssdk.services.batch.model.JobDefinition
import software.amazon.awssdk.services.batch.model.JobDefinitionType
import software.amazon.awssdk.services.batch.model.JobDetail
import software.amazon.awssdk.services.batch.model.JobStatus
import software.amazon.awssdk.services.batch.model.JobTimeout
import software.amazon.awssdk.services.batch.model.KeyValuePair
import software.amazon.awssdk.services.batch.model.LogConfiguration
import software.amazon.awssdk.services.batch.model.MountPoint
import software.amazon.awssdk.services.batch.model.NetworkConfiguration
import software.amazon.awssdk.services.batch.model.PlatformCapability
import software.amazon.awssdk.services.batch.model.RegisterJobDefinitionRequest
import software.amazon.awssdk.services.batch.model.RegisterJobDefinitionResponse
import software.amazon.awssdk.services.batch.model.ResourceRequirement
import software.amazon.awssdk.services.batch.model.ResourceType
import software.amazon.awssdk.services.batch.model.RetryStrategy
import software.amazon.awssdk.services.batch.model.RuntimePlatform
import software.amazon.awssdk.services.batch.model.SubmitJobRequest
import software.amazon.awssdk.services.batch.model.SubmitJobResponse
import software.amazon.awssdk.services.batch.model.TerminateJobRequest
import software.amazon.awssdk.services.batch.model.Volume
/**
 * Implements a task handler for AWS Batch jobs
 */
// note: do not declare this class as `CompileStatic` otherwise the proxy is not get invoked
@Slf4j
class AwsBatchTaskHandler extends TaskHandler implements BatchHandler<String,JobDetail>, FusionAwareTask {

    private final Path exitFile

    private final Path wrapperFile

    private final Path outputFile

    private final Path errorFile

    private final Path logFile

    private final Path scriptFile

    private final Path inputFile

    private final Path traceFile

    private AwsBatchExecutor executor

    private BatchClient client

    private volatile String jobId

    private volatile String taskArn

    private String queueName

    private CloudMachineInfo machineInfo

    private Map<String,String> environment = Map<String,String>.of()

    final static private Map<String,String> jobDefinitions = [:]

    final static private List<String> MISCONFIGURATION_REASONS = List.of(
        "MISCONFIGURATION:JOB_RESOURCE_REQUIREMENT",
        "MISCONFIGURATION:COMPUTE_ENVIRONMENT_MAX_RESOURCE"
    )

    /**
     * Batch context shared between multiple task handlers
     */
    private BatchContext<String,JobDetail> context

    @TestOnly
    protected AwsBatchTaskHandler() {}

    /**
     * Create a new Batch task handler
     *
     * @param task The {@link nextflow.processor.TaskRun} descriptor of the task to run
     * @param executor The {@link AwsBatchExecutor} instance
     */
    AwsBatchTaskHandler(TaskRun task, AwsBatchExecutor executor) {
        super(task)
        this.executor = executor
        this.client = executor.client
        this.environment = SysEnv.get()
        this.logFile = task.workDir.resolve(TaskRun.CMD_LOG)
        this.scriptFile = task.workDir.resolve(TaskRun.CMD_SCRIPT)
        this.inputFile =  task.workDir.resolve(TaskRun.CMD_INFILE)
        this.outputFile = task.workDir.resolve(TaskRun.CMD_OUTFILE)
        this.errorFile = task.workDir.resolve(TaskRun.CMD_ERRFILE)
        this.exitFile = task.workDir.resolve(TaskRun.CMD_EXIT)
        this.wrapperFile = task.workDir.resolve(TaskRun.CMD_RUN)
        this.traceFile = task.workDir.resolve(TaskRun.CMD_TRACE)
    }

    protected String getJobId() { jobId }

    /**
     * @return An instance of {@link AwsOptions} holding Batch specific settings
     */
    protected AwsOptions getAwsOptions() { executor.getAwsOptions() }

    /**
     * Set the batch collector object. This has not to be confused AWSBatch.
     * It is needed to aggregate multiple API requests to a single remote
     * invocation. It can be implemented by any executor, not just AWSBatch.
     *
     * @param context The {@link BatchContext} object to be used
     */
    void batch( BatchContext<String,JobDetail> context ) {
        if( jobId ) {
            context.collect(jobId)
            this.context = context
        }
    }

    private String jobIdsToString(Collection items) {
        final MAX=10
        final sz = items.size()
        items.size()<=MAX ? items.join(', ').toString() : items.take(MAX).join(', ').toString() + ", ... other ${sz-MAX} omitted"
    }

    /**
     * Retrieve Batch job status information
     *
     * @param jobId The Batch job ID
     * @return The associated {@link JobDetail} object or {@code null} if no information is found
     */
    protected JobDetail describeJob(String jobId) {

        Collection batchIds
        if( context ) {
            // check if this response is cached in the batch collector
            if( context.contains(jobId) ) {
                log.trace "[AWS BATCH] hit cache for describe job=$jobId"
                return context.get(jobId)
            }
            log.trace "[AWS BATCH] missed cache for describe job=$jobId"
            // get next 100 job ids for which it's required to check the status
            batchIds = context.getBatchFor(jobId, 100)
        }
        else {
            batchIds = [jobId]
        }

        // retrieve the status for the specified job and along with the next batch
        log.trace "[AWS BATCH] requesting describe jobs=${jobIdsToString(batchIds)}"
        final request = DescribeJobsRequest.builder().jobs(batchIds).build()
        final resp = client.describeJobs(request)
        if( !resp || !resp.jobs() ) {
            log.debug "[AWS BATCH] cannot retrieve running status for job=$jobId"
            return null
        }

        JobDetail result=null
        for( JobDetail entry : resp.jobs() ) {
            // cache the response in the batch collector
            context?.put( entry.jobId(), entry )
            // return the job detail for the specified job
            if( entry.jobId() == jobId )
                result = entry
        }
        if( !result ) {
            log.debug "[AWS BATCH] cannot find running status for job=$jobId"
        }
        else {
            log.trace "[AWS BATCH] Job id=$jobId details=$result"
        }

        return result
    }

    /**
     * {@inheritDoc}
     */
    @Override
    boolean checkIfRunning() {
        if( !jobId || !isSubmitted() )
            return false
        final job = describeJob(jobId)
        final result = job?.status() in [JobStatus.RUNNING, JobStatus.SUCCEEDED, JobStatus.FAILED]
        if( result )
            this.status = TaskStatus.RUNNING
        else
            checkIfUnschedulable(job)
        // fetch the task arn
        if( !taskArn )
            taskArn = job?.container()?.taskArn()
        return result
    }

    protected void checkIfUnschedulable(JobDetail job) {
        if( job ) try {
            checkIfUnschedulable0(job)
        }
        catch (Throwable e) {
            log.warn "Unable to check if job is unschedulable - ${e.message}", e
        }
    }

    private void checkIfUnschedulable0(JobDetail job) {
        final reason = errReason(job)
        if( MISCONFIGURATION_REASONS.any((it) -> reason.contains(it)) ) {
            final msg = "unschedulable AWS Batch job ${jobId} (${task.lazyName()}) - $reason"
            // If indicated in aws.batch config kill the job an produce a failure
            if( executor.awsOptions.terminateUnschedulableJobs() ){
                log.warn("Terminating ${jobId}")
                kill()
                task.error = new ProcessException("Unschedulable AWS Batch job ${jobId} - $reason")
                status = TaskStatus.COMPLETED
            }
            else {
                log.warn "Detected $msg"
            }
        }
    }

    protected String errReason(JobDetail job){
        if(!job)
            return "(unknown)"
        final result = new ArrayList(2)
        if( job.statusReason() )
            result.add(job.statusReason())
        final AttemptContainerDetail container = job.attempts() ? job.attempts()[-1].container() : null
        if( container?.reason() )
            result.add(container.reason())
        return result.join(' - ')
    }

    /**
     * {@inheritDoc}
     */
    @Override
    boolean checkIfCompleted() {
        assert jobId
        if( isCompleted() ) {
            //Task can be marked as completed before running by unschedulable reason. Return true
            return true
        }
        if( !isRunning() )
            return false
        final job = describeJob(jobId)
        final done = job?.status() in [JobStatus.SUCCEEDED, JobStatus.FAILED]
        if( done ) {
            // take the exit code of the container, if 0 (successful) or missing
            // take the exit code from the `.exitcode` file create by nextflow
            // the rationale of this is that, in case of error, the exit code return
            // by the batch API is more reliable.
            task.exitStatus = job.container().exitCode() ?: readExitFile()
            // finalize the task
            task.stdout = outputFile
            if( job?.status() == JobStatus.FAILED || task.exitStatus==Integer.MAX_VALUE ) {
                final reason = errReason(job)
                // retry all CannotPullContainer errors apart when it does not exist or cannot be accessed
                final unrecoverable = reason.contains('CannotPullContainer') && reason.contains('unauthorized')
                task.error = unrecoverable ? new ProcessUnrecoverableException(reason) : new ProcessException(reason)
                task.stderr = executor.getJobOutputStream(jobId) ?: errorFile
            }
            else {
                task.stderr = errorFile
            }
            status = TaskStatus.COMPLETED
            return true
        }
        return false
    }

    private int readExitFile() {
        try {
            exitFile.text as Integer
        }
        catch( Exception e ) {
            log.debug "[AWS BATCH] Cannot read exit status for task: `${task.lazyName()}` | ${e.message}"
            return Integer.MAX_VALUE
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void killTask() {
        assert jobId
        final targetId = normaliseJobId(jobId)
        if( executor.shouldDeleteJob(targetId)) {
            terminateJob(targetId)
        }
    }

    protected String normaliseJobId(String jobId) {
        if( !jobId )
            return null
        return jobId.contains(':')
                ? jobId.split(':')[0]
                : jobId
    }

    protected void terminateJob(String jobId) {
        log.debug "[AWS BATCH] cleanup = terminating job $jobId"
        final req = TerminateJobRequest.builder()
            .jobId(jobId)
            .reason('Job killed by NF')
            .build()
        final batch = bypassProxy(client)
        executor.reaper.submit({
            final resp = batch.terminateJob(req)
            log.debug "[AWS BATCH] cleanup = killing job $jobId"
        })
    }

    @Override
    void prepareLauncher() {
        createTaskWrapper().build()
    }

    /**
     * {@inheritDoc}
     */
    @Override
    void submit() {
        /*
         * create submit request
         */
        final req = newSubmitRequest(task)
        log.trace "[AWS BATCH] new job request > $req"

        /*
         * submit the task execution
         */
        // note use the real client object because this method
        // is supposed to be invoked by the thread pool
        final resp = submit0(bypassProxy(client), req)
        updateStatus(resp.jobId(), req.jobQueue())
        log.debug "[AWS BATCH] Process `${task.lazyName()}` submitted > job=$jobId; work-dir=${task.getWorkDirStr()}"
    }

    /*
     * note: this method cannot be 'private' otherwise subclasses (xpack) will fail invoking it
     */
    protected void updateStatus(String jobId, String queueName) {
        if( task instanceof TaskArrayRun ) {
            // update status for children tasks
            for( int i=0; i<task.children.size(); i++ ) {
                final handler = task.children[i] as AwsBatchTaskHandler
                final arrayTaskId = executor.getArrayTaskId(jobId, i)
                handler.updateStatus(arrayTaskId, queueName)
            }
        }
        else {
            this.jobId = jobId
            this.queueName = queueName
            this.status = TaskStatus.SUBMITTED
        }
    }

    protected BashWrapperBuilder createTaskWrapper() {
        return fusionEnabled()
                ? fusionLauncher()
                : new AwsBatchScriptLauncher(task.toTaskBean(), getAwsOptions())
    }

    protected void buildTaskWrapper() {
        createTaskWrapper().build()
    }

    protected BatchClient bypassProxy(BatchClient batch) {
        batch instanceof AwsBatchProxy ? batch.getClient() : batch
    }

    /**
     * Retrieve the queue name to use to submit the task execution
     *
     * @param task The {@link TaskRun} object to be executed
     * @return The Batch queue name defined by this job execution
     */
    protected String getJobQueue(TaskRun task) {
        final queue = task.config.queue?.toString()
        if( !queue )
            throw new ProcessUnrecoverableException("Missing AWS Batch job queue -- provide it by using the process `queue` directive")

        return queue
    }

    /**
     * Get the Batch job definition name used to run the specified task
     *
     * @param task The {@link TaskRun} object to be executed
     * @return The Batch job definition name defined by this job execution
     */
    protected String getJobDefinition(TaskRun task) {
        final container = task.getContainer()
        if( !container )
            throw new ProcessUnrecoverableException("Invalid AWS Batch job definition -- provide a Docker image name or a Batch job definition name")

        if( container.startsWith('job-definition://')) {
            return container.substring(17)
        }

        return resolveJobDefinition(task)
    }

    /**
     * Maps a docker container image to a Batch job definition name
     *
     * @param container The Docker container image name which need to be used to run the job
     * @return The Batch Job Definition name associated with the specified container
     */
    @CompileStatic
    protected String resolveJobDefinition(TaskRun task) {
        final int DEFAULT_BACK_OFF_BASE = 3
        final int DEFAULT_BACK_OFF_DELAY = 250
        final int MAX_ATTEMPTS = 5
        int attempt=0
        while( true ) {
            try {
                return resolveJobDefinition0(task)
            }
            catch (ClientException e) {
                if( e.statusCode() != 404 || attempt++ > MAX_ATTEMPTS)
                    throw e

                final delay = (Math.pow(DEFAULT_BACK_OFF_BASE, attempt) as long) * DEFAULT_BACK_OFF_DELAY
                log.debug "Got AWS Client exception on Batch resolve job definition - message=$e.message; waiting for ${delay}ms (attempt=$attempt)"
                Thread.sleep(delay)
            }
        }
    }

    @CompileStatic
    protected String resolveJobDefinition0(TaskRun task) {
        final req = makeJobDefRequest(task)
        final container = task.getContainer()
        final token = req.parameters.get('nf-token')
        final jobKey = "$container:$token".toString()
        if( jobDefinitions.containsKey(jobKey) )
            return jobDefinitions[jobKey]

        synchronized(jobDefinitions) {
            if( jobDefinitions.containsKey(jobKey) )
                return jobDefinitions[jobKey]

            def msg
            def name = findJobDef(req.jobDefinitionName, token)
            if( name ) {
                msg = "[AWS BATCH] Found job definition name=$name; container=$container"
            }
            else {
                name = createJobDef(req)
                msg = "[AWS BATCH] Created job definition name=$name; container=$container"
            }
            // log the request
            if( log.isTraceEnabled() )
                log.debug "[AWS BATCH] $msg; request=${req.toString().indent()}"
            else
                log.debug "[AWS BATCH] $msg"

            jobDefinitions[jobKey] = name
            return name
        }
    }

    /**
     * Create a Batch job definition request object for the specified Docker image
     *
     * @param image The Docker container image for which is required to create a Batch job definition
     * @return An instance of {@link RegisterJobDefinitionModel} for the specified Docker image
     */
    @CompileStatic
    protected RegisterJobDefinitionModel makeJobDefRequest(TaskRun task) {
        final uniq = new ArrayList()
        final result = configJobDefRequest(task, uniq)

        // create a job marker uuid
        def hash = computeUniqueToken(uniq)
        result.parameters(['nf-token':hash])

        return result
    }

    protected String computeUniqueToken(List uniq) {
        return CacheHelper.hasher(uniq).hash().toString()
    }

    /**
     * Create and configure the actual RegisterJobDefinitionRequest object
     *
     * @param image
     *      The Docker container image for which is required to create a Batch job definition
     * @param hashingTokens
     *      A list used to collect values that should be used to create a unique job definition Id for the given job request.
     *      It should be used to return such values in the calling context
     * @return
     *      An instance of {@link RegisterJobDefinitionModel} for the specified Docker image
     */
    @CompileStatic
    protected RegisterJobDefinitionModel configJobDefRequest(TaskRun task, List hashingTokens) {
        final image = task.getContainer()
        final name = normalizeJobDefinitionName(image)
        final opts = getAwsOptions()

        final result = new RegisterJobDefinitionModel()
        result.jobDefinitionName(name)
        result.type(JobDefinitionType.CONTAINER)

        // create the container opts based on task config
        final containerOpts = task.getConfig().getContainerOptionsMap()
        final container = createContainerProperties(containerOpts)

        // container definition
        // https://docs.aws.amazon.com/AmazonECS/latest/developerguide/task-cpu-memory-error.html
        final reqCpus = ResourceRequirement.builder().type(ResourceType.VCPU).value('1').build()
        final reqMem = ResourceRequirement.builder().type(ResourceType.MEMORY).value( opts.fargateMode ? '2048' : '1024').build()
        container
                .image(image)
                .command('true')
                // note the actual command, memory and cpus are overridden when the job is executed
                .resourceRequirements( reqCpus, reqMem )

        final jobRole = opts.getJobRole()
        if( jobRole )
            container.jobRoleArn(jobRole)

        if( opts.executionRole )
            container.executionRoleArn(opts.executionRole)
        
        final logsGroup = opts.getLogsGroup()
        if( logsGroup )
            container.logConfiguration(getLogConfiguration(logsGroup, opts.getRegion()))

        if( fusionEnabled() )
            container.privileged(true)

        final mountsMap = new LinkedHashMap( 10)
        final awscli = opts.cliPath
        if( awscli && !opts.fargateMode ) {
            def path = Paths.get(awscli).parent.parent.toString()
            mountsMap.put('aws-cli', "$path:$path:ro")
        }

        int c=0
        final volumes = opts.getVolumes()
        for( String vol : volumes ) {
            mountsMap.put("vol-"+(++c), vol)
        }

        if( mountsMap )
            addVolumeMountsToContainer(mountsMap, container)

        // Fargate specific settings
        if( opts.isFargateMode() ) {
            result.platformCapabilities(List.of(PlatformCapability.FARGATE))
            container.networkConfiguration( NetworkConfiguration.builder().assignPublicIp(AssignPublicIp.ENABLED).build() )
            // use at least 50 GB as disk local storage
            final diskGb = task.config.getDisk()?.toGiga()?.toInteger() ?: 50
            container.ephemeralStorage( EphemeralStorage.builder().sizeInGiB(diskGb).build() )
            // check for arm64 cpu architecture
            if( task.config.getArchitecture()?.arch == 'arm64' )
                container.runtimePlatform(RuntimePlatform.builder().cpuArchitecture('ARM64').build())
        }

        // finally set the container options
        result.containerProperties(container)

        // add to this list all values that has to contribute to the
        // job definition unique name creation
        hashingTokens.add(name)
        hashingTokens.add(container.toString())
        if( containerOpts )
            hashingTokens.add(containerOpts)

        return result
    }

    @Memoized
    LogConfiguration getLogConfiguration(String name, String region) {
        LogConfiguration.builder()
            .logDriver('awslogs')
            .options([
                'awslogs-region': region,
                'awslogs-group': name
            ]).build()
    }

    @CompileStatic
    protected void addVolumeMountsToContainer(Map<String,String> mountsMap, ContainerPropertiesModel container) {
        final mounts = new ArrayList<MountPoint>(mountsMap.size())
        final volumes = new  ArrayList<Volume>(mountsMap.size())
        for( Map.Entry<String,String> entry : mountsMap.entrySet() ) {
            final mountName = entry.key
            final parts = entry.value.tokenize(':')
            final containerPath = parts[0]
            final hostPath = parts.size()>1 ? parts[1] : containerPath
            final readOnly = parts.size()>2 ? parts[2]=='ro' : false
            if( parts.size()>3 )
                throw new IllegalArgumentException("Not a valid volume mount syntax: $entry.value")

            def mount = MountPoint.builder()
                    .sourceVolume(mountName)
                    .containerPath(hostPath)
                    .readOnly(readOnly)
                    .build()
            mounts << mount

            def vol = Volume.builder()
                    .name(mountName)
                    .host(Host.builder().sourcePath(containerPath).build())
                    .build()
            volumes << vol
        }

        if( mountsMap ) {
            container.mountPoints(mounts)
            container.volumes(volumes)
        }
    }

    /**
     * Look for a Batch job definition in ACTIVE status for the given name and NF job definition ID
     *
     * @param name The Batch job definition name
     * @param jobId A unique job definition ID generated by NF
     * @return The fully qualified Batch job definition name eg {@code my-job-definition:3}
     */
    protected String findJobDef(String name, String jobId) {
        log.trace "[AWS BATCH] checking job definition with name=$name; jobid=$jobId"
        final req = DescribeJobDefinitionsRequest.builder()
            .jobDefinitionName(name)
            .build()
        // bypass the proxy because this method is invoked during a
        // job submit request that's already in a separate thread pool request
        // therefore it's protected by a TooManyRequestsException
        final res = describeJobDefinitions0(bypassProxy(this.client), req)
        final jobs = res.jobDefinitions()
        if( jobs.size()==0 )
            return null

        def job = jobs.find { JobDefinition it -> it.status() == 'ACTIVE' && it.parameters()?.'nf-token' == jobId  }
        return job ? "$name:${job.revision()}" : null
    }

    /**
     * Create (aka register) a new Batch job definition
     *
     * @param model A {@link RegisterJobDefinitionRequest} representing the Batch jib definition to create
     * @return The fully qualified Batch job definition name eg {@code my-job-definition:3}
     */
    protected String createJobDef(RegisterJobDefinitionModel model) {
        // add nextflow tags
        model.addTagsEntry('nextflow.io/createdAt', Instant.now().toString())
        model.addTagsEntry('nextflow.io/version', BuildInfo.version)
        // create the job def
        final req = model.toBatchRequest()
        final res = createJobDef0(bypassProxy(client), req) // bypass the client proxy! see #1024
        return "${res.jobDefinitionName()}:${res.revision()}"
    }

    /**
     * Make a name string compliant with the Batch job definition format
     *
     * @param name A job name
     * @return A given name formatted to be used as Job definition name
     */
    protected String normalizeJobDefinitionName(String name) {
        if( !name ) return null
        if( !ContainerNameValidator.isValidImageName(name) ) throw new IllegalArgumentException("Invalid container image name: $name")

        def result = name.replaceAll(/[^a-zA-Z0-9\-_]+/,'-')
        // Batch job definition length cannot exceed 128 characters
        // take first 40 chars + add a unique MD5 hash (32 chars)
        if( result.length()>125 ) {
            final hash = name.md5()
            result = result.substring(0,40) + '-' + hash
        }

        return "nf-" + result
    }

    protected List<String> classicSubmitCli() {
        return executor.getLaunchCommand(task.getWorkDirStr())
    }

    protected List<String> getSubmitCommand() {
        // final launcher command
        return fusionEnabled()
                ? fusionSubmitCli()
                : classicSubmitCli()
    }

    protected int maxSpotAttempts() {
        final result = executor.awsOptions.maxSpotAttempts
        if( result )
            return result
        // when fusion snapshot is enabled max attempt should be > 0
        // to enable to allow snapshot retry the job execution in a new ec2 instance
        return fusionEnabled() && fusionConfig().snapshotsEnabled() ? 5 : 0
    }

    protected String getJobName(TaskRun task) {
        final result = prependWorkflowPrefix(task.name, environment)
        return normalizeJobName(result)
    }

    /**
     * Create a new Batch job request for the given NF {@link TaskRun}
     *
     * @param task A {@link TaskRun} to be executed as Batch job
     * @return A {@link software.amazon.awssdk.services.batch.model.SubmitJobRequest} instance representing the Batch job to submit
     */
    protected SubmitJobRequest newSubmitRequest(TaskRun task) {

        /*
         * create the request object
         */
        final opts = getAwsOptions()
        final labels = task.config.getResourceLabels()
        final builder = SubmitJobRequest.builder()
        builder.jobName(getJobName(task))
        builder.jobQueue(getJobQueue(task))
        builder.jobDefinition(getJobDefinition(task))
        if( labels ) {
            final tags = validateAwsBatchLabels(labels)
            builder.tags(tags)
            builder.propagateTags(true)
        }
        // set the share identifier
        if( opts.shareIdentifier ) {
            builder.shareIdentifier(opts.shareIdentifier)
            builder.schedulingPriorityOverride(opts.schedulingPriority)
        }

        /*
         * retry on spot reclaim
         * https://aws.amazon.com/blogs/compute/introducing-retry-strategies-for-aws-batch/
         */
        final attempts = maxSpotAttempts()
        if( attempts>0 ) {
            // retry the job when an Ec2 instance is terminate
            final cond1 = EvaluateOnExit.builder().action('RETRY').onStatusReason('Host EC2*').build()
            // the exit condition prevent to retry for other reason and delegate
            // instead to nextflow error strategy the handling of the error
            final cond2 = EvaluateOnExit.builder().action('EXIT').onReason('*').build()
            final retry = RetryStrategy.builder()
                    .attempts( attempts )
                    .evaluateOnExit(cond1, cond2)
                    .build()
            builder.retryStrategy(retry)
        }

        // set task timeout
        final time = task.config.getTime()
        if( time ) {
            def secs = time.toSeconds() as Integer
            if( secs < 60 ) {
                secs = 60   // Batch minimal allowed timeout is 60 seconds
            }
            builder.timeout(JobTimeout.builder().attemptDurationSeconds(secs).build())
        }

        // set the actual command
        final resources = new ArrayList<ResourceRequirement>(5)
        final container = ContainerOverrides.builder()
        container.command(getSubmitCommand())
        // set the task memory
        final cpus = task.config.getCpus()
        final mem = task.config.getMemory()
        if( mem ) {
            final mega = opts.fargateMode ? normaliseFargateMem(cpus, mem) : mem.toMega()
            if( mega >= 4 )
                resources << ResourceRequirement.builder().type(ResourceType.MEMORY).value(mega.toString()).build()
            else
                log.warn "Ignoring task ${task.lazyName()} memory directive: ${task.config.getMemory()} -- AWS Batch job memory request cannot be lower than 4 MB"
        }
        // set the task cpus
        if( cpus > 1 )
            resources << ResourceRequirement.builder().type(ResourceType.VCPU).value(task.config.getCpus().toString()).build()

        final accelerator = task.config.getAccelerator()
        if( accelerator ) {
            if( accelerator.type )
                log.warn1 "Ignoring task ${task.lazyName()} accelerator type: ${accelerator.type} -- AWS Batch doesn't support accelerator type in job definition"
            resources << ResourceRequirement.builder().type(ResourceType.GPU).value(accelerator.request.toString()).build()
        }

        if( resources )
            container.resourceRequirements(resources)

        // set the environment
        def vars = getEnvironmentVars()
        if( vars )
            container.environment(vars)

        builder.containerOverrides(container.build())

        // set the array properties
        if( task instanceof TaskArrayRun ) {
            final arraySize = task.getArraySize()

            if( arraySize > 10_000 )
                throw new IllegalArgumentException("Job arrays on AWS Batch may not have more than 10,000 tasks")

            builder.arrayProperties(ArrayProperties.builder().size(arraySize).build())
        }

        return builder.build()
    }

    /**
     * Validate AWS Batch labels for compliance with AWS naming requirements.
     * This method validates resource labels against AWS Batch tag constraints and
     * handles violations based on the nextflow.enable.strict setting:
     * 
     * - When strict mode is disabled (default): logs warnings for invalid tags but allows them through
     * - When strict mode is enabled: throws ProcessUnrecoverableException for invalid tags
     *
     * AWS Batch tag constraints validated:
     * - Keys and values cannot be null
     * - Maximum key length: 128 characters  
     * - Maximum value length: 256 characters
     * - Allowed characters: letters, numbers, spaces, and: _ . : / = + - @
     *
     * @param labels The original resource labels map
     * @return The labels map (unchanged in validation mode)
     * @throws ProcessUnrecoverableException when strict mode is enabled and labels are invalid
     */
    protected Map<String, String> validateAwsBatchLabels(Map<String, String> labels) {
        if (!labels) return labels

        final strictMode = executor.session.config.navigate('nextflow.enable.strict', false)
        final violations = []
        final result = new HashMap<String, String>()

        for (Map.Entry<String, String> entry : labels.entrySet()) {
            final key = entry.getKey()
            final value = entry.getValue()

            // Check for null keys or values and filter them out (not validation violations)
            if (key == null) {
                log.warn "AWS Batch label dropped due to null key: key=null, value=${value}"
                continue
            }
            if (value == null) {
                log.warn "AWS Batch label dropped due to null value: key=${key}, value=null"
                continue
            }

            final keyStr = key.toString()
            final valueStr = value.toString()
            
            // Validate key length
            if (keyStr.length() > 128) {
                violations << "Label key exceeds 128 characters: '${keyStr}' (${keyStr.length()} chars)"
            }
            
            // Validate value length  
            if (valueStr.length() > 256) {
                violations << "Label value exceeds 256 characters: '${keyStr}' = '${valueStr}' (${valueStr.length()} chars)"
            }
            
            // Validate key characters
            if (!isValidAwsBatchTagString(keyStr)) {
                violations << "Label key contains invalid characters: '${keyStr}' - only letters, numbers, spaces, and _ . : / = + - @ are allowed"
            }
            
            // Validate value characters
            if (!isValidAwsBatchTagString(valueStr)) {
                violations << "Label value contains invalid characters: '${keyStr}' = '${valueStr}' - only letters, numbers, spaces, and _ . : / = + - @ are allowed"
            }

            // Add valid entries to result
            result[keyStr] = valueStr
        }

        // Handle violations based on strict mode (but only for constraint violations, not null filtering)
        if (violations) {
            final message = "AWS Batch tag validation failed:\n${violations.collect{ '  - ' + it }.join('\n')}"
            if (strictMode) {
                throw new ProcessUnrecoverableException(message)
            } else {
                log.warn "${message}\nTags will be used as-is but may cause AWS Batch submission failures"
            }
        }

        return result
    }

    /**
     * Check if a string contains only characters allowed in AWS Batch tags.
     * AWS Batch allows: letters, numbers, spaces, and: _ . : / = + - @
     * 
     * @param input The string to validate
     * @return true if the string contains only valid characters
     */
    protected boolean isValidAwsBatchTagString(String input, int maxLength = 128) {
        if (!input) return false
        if (input.length() > maxLength) return false
        return input ==~ /^[a-zA-Z0-9\s_.\:\/=+\-@]*$/
    }

    /**
     * @return The list of environment variables to be defined in the Batch job execution context
     */
    protected List<KeyValuePair> getEnvironmentVars() {
        List<KeyValuePair> vars = []
        if( this.environment?.containsKey('NXF_DEBUG') )
            vars << KeyValuePair.builder().name('NXF_DEBUG').value(this.environment['NXF_DEBUG']).build()
        if( this.getAwsOptions().retryMode && this.getAwsOptions().retryMode in AwsOptions.VALID_RETRY_MODES)
            vars << KeyValuePair.builder().name('AWS_RETRY_MODE').value(this.getAwsOptions().retryMode).build()
        if( this.getAwsOptions().maxTransferAttempts ) {
            vars << KeyValuePair.builder().name('AWS_MAX_ATTEMPTS').value(this.getAwsOptions().maxTransferAttempts as String).build()
            vars << KeyValuePair.builder().name('AWS_METADATA_SERVICE_NUM_ATTEMPTS').value(this.getAwsOptions().maxTransferAttempts as String).build()
        }
        if( fusionEnabled() ) {
            for(Map.Entry<String,String> it : fusionLauncher().fusionEnv()) {
                vars << KeyValuePair.builder().name(it.key).value(it.value).build()
            }
        }
        return vars
    }

    /**
     * @return The launcher script file {@link Path}
     */
    protected Path getWrapperFile() { wrapperFile }

    /**
     * @return The launcher log file {@link Path}
     */
    protected Path getLogFile() { logFile }

    /**
     * Remove invalid characters from a job name string
     *
     * @param name A job name containing possible invalid character
     * @return A job name without invalid characters
     */
    protected String normalizeJobName(String name) {
        def result = name.replaceAll(' ','_').replaceAll(/[^a-zA-Z0-9_-]/,'')
        result.size()>128 ? result.substring(0,128) : result
    }

    protected CloudMachineInfo getMachineInfo() {
        if( machineInfo )
            return machineInfo
        if( queueName && taskArn && executor.awsOptions.fetchInstanceType ) {
            machineInfo = executor.getMachineInfoByQueueAndTaskArn(queueName, taskArn)
            log.trace "[AWS BATCH] jobId=$jobId; queue=$queueName; task=$taskArn => machineInfo=$machineInfo"
        }
        return machineInfo
    }

    TraceRecord getTraceRecord() {
        def result = super.getTraceRecord()
        result.put('native_id', jobId)
        result.machineInfo = getMachineInfo()
        return result
    }

    // -- helpers

    static private SubmitJobResponse submit0(BatchClient client, SubmitJobRequest req) {
        try {
            return client.submitJob(req)
        }
        catch (BatchException e) {
            if( e.awsErrorDetails().sdkHttpResponse().statusCode() >= 500 )
                // raise a process exception so that nextflow can try to recover it
                throw new ProcessSubmitException("Failed to submit job: ${req.jobName()} - Reason: ${e.awsErrorDetails().errorCode()}", e)
            else
                // status code < 500 are not expected to be recoverable, just throw it again
                throw e
        }
    }

    static private DescribeJobDefinitionsResponse describeJobDefinitions0(BatchClient client, DescribeJobDefinitionsRequest req) {
        try {
            client.describeJobDefinitions(req)
        }
        catch (BatchException e) {
            if( e.awsErrorDetails().sdkHttpResponse().statusCode() >= 500 )
                // raise a process exception so that nextflow can try to recover it
                throw new ProcessSubmitException("Failed to describe job definitions: ${req.jobDefinitions()} - Reason: ${e.awsErrorDetails().errorCode()}", e)
            else
                // status code < 500 are not expected to be recoverable, just throw it again
                throw e
        }
    }

    static private RegisterJobDefinitionResponse createJobDef0(BatchClient client, RegisterJobDefinitionRequest req) {
        try {
            return client.registerJobDefinition(req)
        }
        catch (BatchException e) {
            if( e.awsErrorDetails().sdkHttpResponse().statusCode() >= 500 )
                // raise a process exception so that nextflow can try to recover it
                throw new ProcessSubmitException("Failed to register job definition: ${req.jobDefinitionName()} - Reason: ${e.awsErrorDetails().errorCode()}", e)
            else
                // status code < 500 are not expected to be recoverable, just throw it again
                throw e
        }
    }

    @Canonical
    static class MemSlot {
        int min
        int max
        int step

        static ofGiga(int min, int max, int step) {
            new MemSlot(min *1024, max *1024, step *1024)
        }
    }

    static final Map<Integer, MemSlot> FARGATE_MEM = [1 : MemSlot.ofGiga(2,8,1),
                                                      2 : MemSlot.ofGiga(4, 16, 1),
                                                      4 : MemSlot.ofGiga(8, 30, 1),
                                                      8 : MemSlot.ofGiga(16,60, 4),
                                                      16: MemSlot.ofGiga(32, 120, 8) ]

    protected long normaliseFargateMem(Integer cpus, MemoryUnit mem) {
        final mega = mem.toMega()
        final slot = FARGATE_MEM.get(cpus)
        if( slot==null )
            throw new ProcessUnrecoverableException("Requirement of $cpus CPUs is not allowed by Fargate -- Check process with name '${task.lazyName()}'")
        if( mega <slot.min ) {
            log.warn "Process '${task.lazyName()}' memory requirement of ${mem} is below the minimum allowed by Fargate of ${MemoryUnit.of(mega+'MB')}"
            return slot.min
        }
        if( mega >slot.max ) {
            log.warn "Process '${task.lazyName()}' memory requirement of ${mem} is above the maximum allowed by Fargate of ${MemoryUnit.of(mega+'MB')}"
            return slot.max
        }
        return ceilDiv(mega, slot.step) * slot.step
    }

    static private long ceilDiv(long x, long y){
        return -Math.floorDiv(-x,y);
    }
}

