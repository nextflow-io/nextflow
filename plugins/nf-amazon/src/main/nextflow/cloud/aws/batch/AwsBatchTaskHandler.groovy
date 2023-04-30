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
import com.amazonaws.services.batch.model.AWSBatchException
import com.amazonaws.services.batch.model.AttemptContainerDetail
import com.amazonaws.services.batch.model.DescribeJobsRequest
import com.amazonaws.services.batch.model.DescribeJobsResult
import com.amazonaws.services.batch.model.JobDetail
import com.amazonaws.services.batch.model.SubmitJobRequest
import com.amazonaws.services.batch.model.SubmitJobResult
import com.amazonaws.services.batch.model.TerminateJobRequest
import groovy.transform.CompileDynamic
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.cloud.types.CloudMachineInfo
import nextflow.exception.ProcessSubmitException
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

    @Override
    AWSBatch getClient() { executor.client }

    @Override
    AwsOptions getAwsOptions() { executor.getAwsOptions() }

    protected Path getWrapperFile() { wrapperFile }

    protected Path getLogFile() { logFile }

    protected String getJobId() { jobId }

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

    @Override
    Path prepareLauncher() {
        createTaskWrapper().build()
    }

    protected BashWrapperBuilder createTaskWrapper() {
        return fusionEnabled()
                ? fusionLauncher()
                : new AwsBatchScriptLauncher(task.toTaskBean(), getAwsOptions())
    }

    @Override
    List<String> getSubmitCommand() {
        return fusionEnabled()
                ? fusionSubmitCli()
                : classicSubmitCli()
    }

    protected List<String> classicSubmitCli() {
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
    void submit() {
        if( arraySubmitter ) {
            arraySubmitter.collect(this)
            return
        }

        // -- create the submit request
        final request = newSubmitRequest(task)
        log.trace "[AWS BATCH] new job request > $request"

        // -- submit the job
        // use the real client object because this method
        // should be invoked by the thread pool
        final response = submitJobRequest0(bypassProxy(client), request)
        this.jobId = response.jobId
        this.queueName = request.getJobQueue()
        this.status = TaskStatus.SUBMITTED
        log.debug "[AWS BATCH] submitted > job=$jobId; work-dir=${task.getWorkDirStr()}"
    }

    static private SubmitJobResult submitJobRequest0(AWSBatch client, SubmitJobRequest request) {
        try {
            return client.submitJob(request)
        }
        catch( AWSBatchException e ) {
            if( e.statusCode >= 500 )
                // raise a process exception so that nextflow can try to recover it
                throw new ProcessSubmitException("Failed to submit job: ${request.jobName} - Reason: ${e.errorCode}", e)
            else
                // status code < 500 are not expected to be recoverable, just throw it again
                throw e
        }
    }

    void setJobId(String jobId) {
        this.jobId = jobId
    }

    void setQueueName(String queueName) {
        this.queueName = queueName
    }

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

    private String jobIdsToString(Collection items) {
        final MAX = 10
        final sz = items.size()
        items.size() <= MAX
            ? items.join(', ').toString()
            : items.take(MAX).join(', ').toString() + ", ... other ${sz - MAX} omitted"
    }

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
    TraceRecord getTraceRecord() {
        def record = super.getTraceRecord()
        record.put('native_id', jobId)
        record.machineInfo = getMachineInfo()
        return record
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

}

