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

import java.nio.file.Path
import java.util.concurrent.TimeUnit

import com.amazonaws.services.batch.AWSBatch
import com.amazonaws.services.batch.model.AWSBatchException
import com.amazonaws.services.ecs.model.AccessDeniedException
import com.amazonaws.services.logs.model.ResourceNotFoundException
import groovy.transform.CompileDynamic
import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.cloud.aws.AwsClientFactory
import nextflow.cloud.aws.config.AwsConfig
import nextflow.cloud.aws.nio.S3Path
import nextflow.cloud.types.CloudMachineInfo
import nextflow.exception.AbortOperationException
import nextflow.executor.Executor
import nextflow.executor.TaskArrayExecutor
import nextflow.extension.FilesEx
import nextflow.fusion.FusionHelper
import nextflow.processor.ParallelPollingMonitor
import nextflow.processor.TaskHandler
import nextflow.processor.TaskMonitor
import nextflow.processor.TaskRun
import nextflow.util.Duration
import nextflow.util.Escape
import nextflow.util.RateUnit
import nextflow.util.ServiceName
import nextflow.util.ThreadPoolHelper
import nextflow.util.ThrottlingExecutor
import org.pf4j.ExtensionPoint
/**
 * AWS Batch executor
 * https://aws.amazon.com/batch/
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@ServiceName('awsbatch')
@CompileStatic
class AwsBatchExecutor extends Executor implements ExtensionPoint, TaskArrayExecutor {

    /**
     * Proxy to throttle AWS batch client requests
     */
    @PackageScope
    private AwsBatchProxy client

    /** Helper class to resolve Batch related metadata */
    private AwsBatchHelper helper

    /**
     * executor service to throttle service requests
     */
    private ThrottlingExecutor submitter

    /**
     * Executor service to throttle cancel requests
     */
    private ThrottlingExecutor reaper

    /**
     * A S3 path where executable scripts need to be uploaded
     */
    private Path remoteBinDir = null

    private AwsOptions awsOptions

    private final Set<String> deletedJobs = new HashSet<>()

    AwsOptions getAwsOptions() {  awsOptions  }

    /**
     * @return {@code true} to signal containers are managed directly the AWS Batch service
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
        if( !(workDir instanceof S3Path) ) {
            session.abort()
            throw new AbortOperationException("When using `$name` executor an S3 bucket must be provided as working directory using either the `-bucket-dir` or `-work-dir` command line option")
        }
    }

    protected void validatePathDir() {
        def path = session.config.navigate('env.PATH')
        if( path ) {
            log.warn "Environment PATH defined in config file is ignored by AWS Batch executor"
        }
    }

    protected void uploadBinDir() {
        /*
         * upload local binaries
         */
        if( session.binDir && !session.binDir.empty() && !session.disableRemoteBinDir ) {
            def s3 = getTempDir()
            log.info "Uploading local `bin` scripts folder to ${s3.toUriString()}/bin"
            remoteBinDir = FilesEx.copyTo(session.binDir, s3)
        }
    }

    protected void createAwsClient() {
        /*
         * retrieve config and credentials and create AWS client
         */
        final driver = new AwsClientFactory(new AwsConfig(session.config.aws as Map))

        /*
         * create a proxy for the aws batch client that manages the request throttling
         */
        client = new AwsBatchProxy(driver.getBatchClient(), submitter)
        helper = new AwsBatchHelper(client, driver)
        // create the options object
        awsOptions = new AwsOptions(this)
        log.debug "[AWS BATCH] Executor ${awsOptions.fargateMode ? '(FARGATE mode) ' : ''}options=$awsOptions"
    }

    /**
     * Initialise the AWS batch executor.
     */
    @Override
    protected void register() {
        super.register()
        validateWorkDir()
        validatePathDir()
        uploadBinDir()
        createAwsClient()
    }

    @PackageScope
    Path getRemoteBinDir() {
        remoteBinDir
    }

    @PackageScope
    AWSBatch getClient() {
        client
    }

    /**
     * @return The monitor instance that handles AWS batch tasks
     */
    @Override
    protected TaskMonitor createTaskMonitor() {

        // create the throttling executor
        // note this is invoke only the very first time a AWS Batch executor is created
        // therefore it's safe to assign to a static attribute
        submitter = createExecutorService('AWSBatch-executor')

        reaper = createExecutorService('AWSBatch-reaper')

        final pollInterval = session.getPollInterval(name, Duration.of('10 sec'))
        final dumpInterval = session.getMonitorDumpInterval(name)
        final capacity = session.getQueueSize(name, 1000)

        final def params = [
                name: name,
                session: session,
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
     * @return A {@link AwsBatchTaskHandler} for the given task
     */
    @Override
    TaskHandler createTaskHandler(TaskRun task) {
        assert task
        assert task.workDir
        log.trace "[AWS BATCH] Launching process > ${task.name} -- work folder: ${task.workDirStr}"
        new AwsBatchTaskHandler(task, this)
    }

    private static final List<Integer> RETRYABLE_STATUS = [429, 500, 502, 503, 504]

    /**
     * @return Creates a {@link ThrottlingExecutor} service to throttle
     * the API requests to the AWS Batch service.
     */
    private ThrottlingExecutor createExecutorService(String name) {

        // queue size can be overridden by submitter options below
        final qs = 5_000
        final limit = session.getExecConfigProp(name,'submitRateLimit','50/s') as String
        final size = Runtime.runtime.availableProcessors() * 5

        final opts = new ThrottlingExecutor.Options()
                .retryOn { Throwable t -> t instanceof AWSBatchException && (t.errorCode=='TooManyRequestsException' || t.statusCode in RETRYABLE_STATUS) }
                .onFailure { Throwable t -> session?.abort(t) }
                .onRateLimitChange { RateUnit rate -> logRateLimitChange(rate) }
                .withRateLimit(limit)
                .withQueueSize(qs)
                .withPoolSize(size)
                .withKeepAlive(Duration.of('1 min'))
                .withAutoThrottle(true)
                .withMaxRetries(10)
                .withOptions( getConfigOpts() )
                .withPoolName(name)

        ThrottlingExecutor.create(opts)
    }

    @CompileDynamic
    protected Map getConfigOpts() {
        session.config?.executor?.submitter as Map
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

    boolean shouldDeleteJob(String jobId) {
        if( jobId in deletedJobs ) {
            // if the job is already in the list if has been already deleted
            return false
        }
        synchronized (deletedJobs) {
            // add the job id to the set of deleted jobs, if it's a new id, the `add` method
            // returns true therefore the job should be deleted
            return deletedJobs.add(jobId)
        }
    }

    CloudMachineInfo getMachineInfoByQueueAndTaskArn(String queue, String taskArn) {
        try {
            return helper?.getCloudInfoByQueueAndTaskArn(queue, taskArn)
        }
        catch ( AccessDeniedException e ) {
            log.warn "Unable to retrieve AWS Batch instance type | ${e.message}"
            // disable it since user has not permission to access this info
            awsOptions.fetchInstanceType = false
            return null
        }
        catch( Exception e ) {
            log.warn "Unable to retrieve AWS batch instance type for queue=$queue; task=$taskArn | ${e.message}", e
            return null
        }
    }

    String getJobOutputStream(String jobId) {
        try {
            return helper.getTaskLogStream(jobId, awsOptions.getLogsGroup())
        }
        catch (ResourceNotFoundException e) {
            log.debug "Unable to find AWS Cloudwatch logs for Batch Job id=$jobId - ${e.message}"
        }
        catch (Exception e) {
            log.debug "Unable to retrieve AWS Cloudwatch logs for Batch Job id=$jobId | ${e.message}", e
        }
        return null
    }

    @Override
    void shutdown() {
        def tasks = submitter.shutdownNow()
        if( tasks ) log.warn "Execution interrupted -- cleaning up execution pool"
        submitter.awaitTermination(5, TimeUnit.MINUTES)
        // -- finally delete cleanup executor
        // start shutdown process
        reaper.shutdown()
        final waitMsg = "[AWS BATCH] Waiting jobs reaper to complete (%d jobs to be terminated)"
        final exitMsg = "[AWS BATCH] Exiting before jobs reaper thread pool complete -- Some jobs may not be terminated"
        ThreadPoolHelper.await(reaper, Duration.of('60min'), waitMsg, exitMsg)
    }

    @Override
    String getArrayIndexName() { 'AWS_BATCH_JOB_ARRAY_INDEX' }

    @Override
    int getArrayIndexStart() { 0 }

    @Override
    String getArrayTaskId(String jobId, int index) {
        return "${jobId}:${index}"
    }

    @Override
    String getArrayLaunchCommand(String taskDir) {
        if( isFusionEnabled() || isWorkDirDefaultFS() )
            return TaskArrayExecutor.super.getArrayLaunchCommand(taskDir)
        else
            return Escape.cli(getLaunchCommand(taskDir) as String[])
    }

    List<String> getLaunchCommand(String s3WorkDir) {
        // the cmd list to launch it
        final opts = getAwsOptions()
        final cmd = opts.s5cmdPath
            ? s5Cmd(s3WorkDir, opts)
            : s3Cmd(s3WorkDir, opts)
        return ['bash','-o','pipefail','-c', cmd.toString()]
    }

    static String s3Cmd(String workDir, AwsOptions opts) {
        final cli = opts.getAwsCli()
        final debug = opts.debug ? ' --debug' : ''
        final sse = opts.storageEncryption ? " --sse $opts.storageEncryption" : ''
        final kms = opts.storageKmsKeyId ? " --sse-kms-key-id $opts.storageKmsKeyId" : ''
        final requesterPays = opts.requesterPays ? ' --request-payer requester' : ''
        final aws = "$cli s3 cp --only-show-errors${sse}${kms}${debug}${requesterPays}"
        final cmd = "trap \"{ ret=\$?; $aws ${TaskRun.CMD_LOG} ${workDir}/${TaskRun.CMD_LOG}||true; exit \$ret; }\" EXIT; $aws ${workDir}/${TaskRun.CMD_RUN} - | bash 2>&1 | tee ${TaskRun.CMD_LOG}"
        return cmd
    }

    static String s5Cmd(String workDir, AwsOptions opts) {
        final cli = opts.getS5cmdPath()
        final sse = opts.storageEncryption ? " --sse $opts.storageEncryption" : ''
        final kms = opts.storageKmsKeyId ? " --sse-kms-key-id $opts.storageKmsKeyId" : ''
        final requesterPays = opts.requesterPays ? ' --request-payer requester' : ''
        final cmd = "trap \"{ ret=\$?; $cli cp${sse}${kms}${requesterPays} ${TaskRun.CMD_LOG} ${workDir}/${TaskRun.CMD_LOG}||true; exit \$ret; }\" EXIT; $cli cat ${workDir}/${TaskRun.CMD_RUN} | bash 2>&1 | tee ${TaskRun.CMD_LOG}"
        return cmd
    }

}
