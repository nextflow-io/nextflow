/*
 * Copyright 2013-2023, Seqera Labs
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

import com.amazonaws.services.batch.AWSBatch
import com.amazonaws.services.batch.model.AttemptContainerDetail
import com.amazonaws.services.batch.model.DescribeJobsRequest
import com.amazonaws.services.batch.model.DescribeJobsResult
import com.amazonaws.services.batch.model.JobDetail
import com.amazonaws.services.batch.model.TerminateJobRequest
import groovy.transform.CompileDynamic
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.cloud.types.CloudMachineInfo
import nextflow.exception.ProcessUnrecoverableException
import nextflow.executor.BashWrapperBuilder
import nextflow.processor.BatchContext
import nextflow.processor.BatchHandler
import nextflow.processor.TaskHandler
import nextflow.processor.TaskRun
import nextflow.processor.TaskStatus
import nextflow.trace.TraceRecord
/**
 * Implements a task handler for AWS Batch jobs.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class AwsBatchTaskHandler extends TaskHandler implements BatchHandler<String,JobDetail>, SubmitJobAware {

    private final Path exitFile

    private final Path wrapperFile

    private final Path outputFile

    private final Path errorFile

    private final Path logFile

    private final Path scriptFile

    private final Path inputFile

    private final Path traceFile

    private AwsBatchExecutor executor

    private volatile String jobId

    private volatile String taskArn

    private volatile String queueName

    private volatile CloudMachineInfo machineInfo

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
        this.executor = executor
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

    @Override
    AWSBatch getClient() { executor.client }

    /**
     * @return An instance of {@link AwsOptions} holding Batch specific settings
     */
    @Override
    AwsOptions getAwsOptions() { executor.getAwsOptions() }

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
    @CompileDynamic
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
            taskArn = job?.getContainer()?.getTaskArn()
        return result
    }

    protected String errReason(JobDetail job){
        if(!job)
            return "(unknown)"
        final result = new ArrayList(2)
        if( job.statusReason )
            result.add(job.statusReason)
        final AttemptContainerDetail container = job.attempts ? job.attempts[-1].container : null
        if( container?.reason )
            result.add(container.reason)
        return result.join(' - ')
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
            // take the exit code of the container, if 0 (successful) or missing
            // take the exit code from the `.exitcode` file create by nextflow
            // the rationale of this is that, in case of error, the exit code return
            // by the batch API is more reliable.
            task.exitStatus = job.container.exitCode ?: readExitFile()
            // finalize the task
            task.stdout = outputFile
            if( job?.status == 'FAILED' ) {
                task.error = new ProcessUnrecoverableException(errReason(job))
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

    @Override
    void prepareLauncher() {
        createTaskWrapper().build()
    }

    /**
     * {@inheritDoc}
     */
    @Override
    void submit() {
        if( arraySubmitter ) {
            arraySubmitter.collect(this)
            return
        }

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
        final resp = submitJobRequest(bypassProxy(client), req)
        this.jobId = resp.jobId
        this.status = TaskStatus.SUBMITTED
        this.queueName = req.getJobQueue()
        log.debug "[AWS BATCH] submitted > job=$jobId; work-dir=${task.getWorkDirStr()}"
    }

    void setJobId(String jobId) {
        this.jobId = jobId
    }

    void setQueueName(String queueName) {
        this.queueName = queueName
    }

    protected BashWrapperBuilder createTaskWrapper() {
        return fusionEnabled()
                ? fusionLauncher()
                : new AwsBatchScriptLauncher(task.toTaskBean(), getAwsOptions())
    }

    protected List<String> classicSubmitCli() {
        // the cmd list to launch it
        final opts = getAwsOptions()
        final cli = opts.getAwsCli()
        final debug = opts.debug ? ' --debug' : ''
        final sse = opts.storageEncryption ? " --sse $opts.storageEncryption" : ''
        final kms = opts.storageKmsKeyId ? " --sse-kms-key-id $opts.storageKmsKeyId" : ''
        final aws = "$cli s3 cp --only-show-errors${sse}${kms}${debug}"
        final cmd = "trap \"{ ret=\$?; $aws ${TaskRun.CMD_LOG} s3:/${getLogFile()}||true; exit \$ret; }\" EXIT; $aws s3:/${getWrapperFile()} - | bash 2>&1 | tee ${TaskRun.CMD_LOG}"
        return ['bash','-o','pipefail','-c', cmd.toString()]
    }

    @Override
    List<String> getSubmitCommand() {
        // final launcher command
        return fusionEnabled()
                ? fusionSubmitCli()
                : classicSubmitCli()
    }

    /**
     * @return The launcher script file {@link Path}
     */
    protected Path getWrapperFile() { wrapperFile }

    /**
     * @return The launcher log file {@link Path}
     */
    protected Path getLogFile() { logFile }

    protected CloudMachineInfo getMachineInfo() {
        if( machineInfo )
            return machineInfo
        if( queueName && taskArn && executor.awsOptions.fetchInstanceType ) {
            machineInfo = executor.getMachineInfoByQueueAndTaskArn(queueName, taskArn)
            log.trace "[AWS BATCH] jobId=$jobId; queue=$queueName; task=$taskArn => machineInfo=$machineInfo"
        }
        return machineInfo
    }

    @Override
    TraceRecord getTraceRecord() {
        def result = super.getTraceRecord()
        result.put('native_id', jobId)
        result.machineInfo = getMachineInfo()
        return result
    }

}

