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

import java.nio.file.Path
import java.nio.file.Paths

import com.amazonaws.services.batch.AWSBatch
import com.amazonaws.services.batch.model.ContainerOverrides
import com.amazonaws.services.batch.model.ContainerProperties
import com.amazonaws.services.batch.model.DescribeJobDefinitionsRequest
import com.amazonaws.services.batch.model.DescribeJobsRequest
import com.amazonaws.services.batch.model.DescribeJobsResult
import com.amazonaws.services.batch.model.Host
import com.amazonaws.services.batch.model.JobDefinition
import com.amazonaws.services.batch.model.JobDefinitionType
import com.amazonaws.services.batch.model.JobDetail
import com.amazonaws.services.batch.model.JobTimeout
import com.amazonaws.services.batch.model.KeyValuePair
import com.amazonaws.services.batch.model.MountPoint
import com.amazonaws.services.batch.model.RegisterJobDefinitionRequest
import com.amazonaws.services.batch.model.ResourceRequirement
import com.amazonaws.services.batch.model.RetryStrategy
import com.amazonaws.services.batch.model.SubmitJobRequest
import com.amazonaws.services.batch.model.TerminateJobRequest
import com.amazonaws.services.batch.model.Volume
import groovy.util.logging.Slf4j
import nextflow.cloud.types.CloudMachineInfo
import nextflow.container.ContainerNameValidator
import nextflow.exception.ProcessUnrecoverableException
import nextflow.executor.BashWrapperBuilder
import nextflow.executor.res.AcceleratorResource
import nextflow.processor.BatchContext
import nextflow.processor.BatchHandler
import nextflow.processor.ErrorStrategy
import nextflow.processor.TaskBean
import nextflow.processor.TaskHandler
import nextflow.processor.TaskRun
import nextflow.processor.TaskStatus
import nextflow.trace.TraceRecord
import nextflow.util.CacheHelper
/**
 * Implements a task handler for AWS Batch jobs
 */
// note: do not declare this class as `CompileStatic` otherwise the proxy is not get invoked
@Slf4j
class AwsBatchTaskHandler extends TaskHandler implements BatchHandler<String,JobDetail> {

    private final Path exitFile

    private final Path wrapperFile

    private final Path outputFile

    private final Path errorFile

    private final Path logFile

    private final Path scriptFile

    private final Path inputFile

    private final Path traceFile

    private TaskBean bean

    private AwsBatchExecutor executor

    private AWSBatch client

    private volatile String jobId

    private volatile String taskArn

    private String queueName

    private CloudMachineInfo machineInfo

    private Map<String,String> environment

    final static private Map<String,String> jobDefinitions = [:]

    /**
     * Batch context shared between multiple task handlers
     */
    private BatchContext<String,JobDetail> context

    /** only for testing purpose -- do not use */
    protected AwsBatchTaskHandler() {}

    /**
     * Create a new Batch task handler
     *
     * @param task The {@link nextflow.processor.TaskRun} descriptor of the task to run
     * @param executor The {@link AwsBatchExecutor} instance
     */
    AwsBatchTaskHandler(TaskRun task, AwsBatchExecutor executor) {
        super(task)
        this.bean = new TaskBean(task)
        this.executor = executor
        this.client = executor.client
        this.environment = System.getenv()

        this.logFile = task.workDir.resolve(TaskRun.CMD_LOG)
        this.scriptFile = task.workDir.resolve(TaskRun.CMD_SCRIPT)
        this.inputFile =  task.workDir.resolve(TaskRun.CMD_INFILE)
        this.outputFile = task.workDir.resolve(TaskRun.CMD_OUTFILE)
        this.errorFile = task.workDir.resolve(TaskRun.CMD_ERRFILE)
        this.exitFile = task.workDir.resolve(TaskRun.CMD_EXIT)
        this.wrapperFile = task.workDir.resolve(TaskRun.CMD_RUN)
        this.traceFile = task.workDir.resolve(TaskRun.CMD_TRACE)
    }


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
        DescribeJobsResult resp = client.describeJobs(new DescribeJobsRequest().withJobs(batchIds))
        if( !resp.getJobs() ) {
            log.debug "[AWS BATCH] cannot retrieve running status for job=$jobId"
            return null
        }

        JobDetail result=null
        for( JobDetail entry : resp.jobs ) {
            // cache the response in the batch collector
            context?.put( entry.jobId, entry )
            // return the job detail for the specified job
            if( entry.jobId == jobId )
                result = entry
        }
        if( !result ) {
            log.debug "[AWS BATCH] cannot find running status for job=$jobId"
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
        final result = job?.status in ['RUNNING', 'SUCCEEDED', 'FAILED']
        if( result )
            this.status = TaskStatus.RUNNING
        // fetch the task arn
        if( !taskArn )
            taskArn = job.getContainer().getTaskArn()
        return result
    }

    /**
     * {@inheritDoc}
     */
    @Override
    boolean checkIfCompleted() {
        assert jobId
        if( !isRunning() )
            return false
        final job = describeJob(jobId)
        final done = job?.status in ['SUCCEEDED', 'FAILED']
        if( done ) {
            // finalize the task
            task.exitStatus = readExitFile()
            task.stdout = outputFile
            task.stderr = errorFile
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
            log.debug "[AWS BATCH] Cannot read exitstatus for task: `$task.name` | ${e.message}"
            return Integer.MAX_VALUE
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    void kill() {
        assert jobId
        log.trace "[AWS BATCH] killing job=$jobId"
        final req = new TerminateJobRequest().withJobId(jobId).withReason('Job killed by NF')
        terminateJob(req)
    }

    protected void terminateJob(TerminateJobRequest req) {

        final batch = bypassProxy(client)
        executor.reaper.submit({
            final resp = batch.terminateJob(req)
            log.debug "[AWS BATCH] killing job=$jobId; response=$resp"
        })
    }

    /**
     * {@inheritDoc}
     */
    @Override
    void submit() {
        /*
         * create task wrapper
         */
        buildTaskWrapper()

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
        final resp = bypassProxy(client).submitJob(req)
        this.jobId = resp.jobId
        this.status = TaskStatus.SUBMITTED
        this.queueName = req.getJobQueue()
        log.debug "[AWS BATCH] submitted > job=$jobId; work-dir=${task.getWorkDirStr()}"
    }

    protected BashWrapperBuilder createTaskWrapper() {
        new AwsBatchScriptLauncher(bean,getAwsOptions())
    }

    protected void buildTaskWrapper() {
        createTaskWrapper().build()
    }

    protected AWSBatch bypassProxy(AWSBatch batch) {
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

        resolveJobDefinition(container)
    }

    /**
     * Maps a docker container image to a Batch job definition name
     *
     * @param container The Docker container image name which need to be used to run the job
     * @return The Batch Job Definition name associated with the specified container
     */
    protected String resolveJobDefinition(String container) {
        if( jobDefinitions.containsKey(container) )
            return jobDefinitions[container]

        // this could race a race condition on the creation of a job definition
        // with another NF instance using the same container name
        synchronized(jobDefinitions) {
            if( jobDefinitions.containsKey(container) )
                return jobDefinitions[container]

            def req = makeJobDefRequest(container)
            def name = findJobDef(req.jobDefinitionName, req.parameters?.'nf-token')
            if( name ) {
                log.debug "[AWS BATCH] Found job definition name=$name; container=$container"
            }
            else {
                name = createJobDef(req)
                log.debug "[AWS BATCH] Created job definition name=$name; container=$container"
            }

            jobDefinitions[container] = name
            return name
        }
    }

    /**
     * Create a Batch job definition request object for the specified Docker image
     *
     * @param image The Docker container image for which is required to create a Batch job definition
     * @return An instance of {@link com.amazonaws.services.batch.model.RegisterJobDefinitionRequest} for the specified Docker image
     */
    protected RegisterJobDefinitionRequest makeJobDefRequest(String image) {
        final name = normalizeJobDefinitionName(image)
        final result = new RegisterJobDefinitionRequest()
        final opts = getAwsOptions()
        result.setJobDefinitionName(name)
        result.setType(JobDefinitionType.Container)


        // container definition
        final container = new ContainerProperties()
                .withImage(image)
        // note the actual command, memory and cpus are overridden when the job is executed
                .withCommand('true')
                .withMemory(1024)
                .withVcpus(1)

        final jobRole = opts.getJobRole()
        if( jobRole )
            container.setJobRoleArn(jobRole)
        
        final mountsMap = new LinkedHashMap( 10)
        final awscli = opts.cliPath
        if( awscli ) {
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

        // finally set the container options
        result.setContainerProperties(container)

        // create a job marker uuid
        def uuid = CacheHelper.hasher([name, image, awscli, volumes, jobRole]).hash().toString()
        result.setParameters(['nf-token':uuid])

        return result
    }

    protected void addVolumeMountsToContainer(Map<String,String> mountsMap, ContainerProperties container) {
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

            def mount = new MountPoint()
                    .withSourceVolume(mountName)
                    .withContainerPath(hostPath)
                    .withReadOnly(readOnly)
            mounts << mount

            def vol = new Volume()
                    .withName(mountName)
                    .withHost(new Host()
                    .withSourcePath(containerPath))
            volumes << vol
        }

        if( mountsMap ) {
            container.setMountPoints(mounts)
            container.setVolumes(volumes)
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
        final req = new DescribeJobDefinitionsRequest().withJobDefinitionName(name)
        // bypass the proxy because this method is invoked during a
        // job submit request that's already in a separate thread pool request
        // therefore it's protected by a TooManyRequestsException
        final res = bypassProxy(this.client).describeJobDefinitions(req)
        final jobs = res.getJobDefinitions()
        if( jobs.size()==0 )
            return null

        def job = jobs.find { JobDefinition it -> it.status == 'ACTIVE' && it.parameters?.'nf-token' == jobId  }
        return job ? "$name:$job.revision" : null
    }

    /**
     * Create (aka register) a new Batch job definition
     *
     * @param req A {@link RegisterJobDefinitionRequest} representing the Batch jib definition to create
     * @return The fully qualified Batch job definition name eg {@code my-job-definition:3}
     */
    protected String createJobDef(RegisterJobDefinitionRequest req) {
        final res = bypassProxy(client).registerJobDefinition(req) // bypass the client proxy! see #1024
        return "${res.jobDefinitionName}:$res.revision"
    }

    /**
     * Make a name string complaint with the Batch job definition format
     *
     * @param name A job name
     * @return A given name formatted to be used as Job definition name
     */
    protected String normalizeJobDefinitionName(String name) {
        if( !name ) return null
        if( !ContainerNameValidator.isValidImageName(name) ) throw new IllegalArgumentException("Invalid container image name: $name")

        def result = name.replaceAll(/[^a-zA-Z0-9\-_]+/,'-')
        return "nf-" + result
    }

    protected List<String> getSubmitCommand() {
        // the cmd list to launch it
        def opts = getAwsOptions()
        def aws = opts.getAwsCli()
        def cmd = "trap \"{ ret=\$?; $aws s3 cp --only-show-errors ${TaskRun.CMD_LOG} s3:/${getLogFile()}||true; exit \$ret; }\" EXIT; $aws s3 cp --only-show-errors s3:/${getWrapperFile()} - | bash 2>&1 | tee ${TaskRun.CMD_LOG}"
        // final launcher command
        return ['bash','-o','pipefail','-c', cmd.toString() ]
    }

    /**
     * Create a new Batch job request for the given NF {@link TaskRun}
     *
     * @param task A {@link TaskRun} to be executed as Batch job
     * @return A {@link com.amazonaws.services.batch.model.SubmitJobRequest} instance representing the Batch job to submit
     */
    protected SubmitJobRequest newSubmitRequest(TaskRun task) {

        /*
         * create the request object
         */
        def result = new SubmitJobRequest()
        result.setJobName(normalizeJobName(task.name))
        result.setJobQueue(getJobQueue(task))
        result.setJobDefinition(getJobDefinition(task))

        // NF uses `maxRetries` *only* if `retry` error strategy is specified
        // otherwise delegates the the retry to AWS Batch
        if( task.config.getMaxRetries() && task.config.getErrorStrategy() != ErrorStrategy.RETRY ) {
            def retry = new RetryStrategy().withAttempts( task.config.getMaxRetries()+1 )
            result.setRetryStrategy(retry)
        }

        // set task timeout
        final time = task.config.getTime()
        if( time ) {
            def secs = time.toSeconds() as Integer
            if( secs < 60 ) {
                secs = 60   // Batch minimal allowed timeout is 60 seconds
            }
            result.setTimeout(new JobTimeout().withAttemptDurationSeconds(secs))
        }

        // set the actual command
        def container = new ContainerOverrides()
        container.command = getSubmitCommand()
        // set the task memory
        if( task.config.getMemory() )
            container.memory = (int)task.config.getMemory().toMega()
        // set the task cpus
        if( task.config.getCpus() > 1 )
            container.vcpus = task.config.getCpus()

        if( task.config.getAccelerator() ) {
            container.setResourceRequirements(createResList(task.config.getAccelerator()))
        }

        // set the environment
        def vars = getEnvironmentVars()
        if( vars )
            container.setEnvironment(vars)

        result.setContainerOverrides(container)

        return result
    }

    protected List<ResourceRequirement> createResList(AcceleratorResource acc) {
        final res = new ResourceRequirement()
        final type = acc.type ?: 'GPU'
        final count = acc.request?.toString() ?: '1'
        res.setType(type.toUpperCase())
        res.setValue(count)
        final result = new ArrayList(1)
        result.add(res)
        return result
    }

    /**
     * @return The list of environment variables to be defined in the Batch job execution context
     */
    protected List<KeyValuePair> getEnvironmentVars() {
        def vars = []
        if( this.environment?.containsKey('NXF_DEBUG') )
            vars << new KeyValuePair().withName('NXF_DEBUG').withValue(this.environment['NXF_DEBUG'])
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
        def result = name.replaceAll(' ','_').replaceAll(/[^a-zA-Z0-9_]/,'')
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
}

