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

import com.amazonaws.services.batch.AWSBatch
import com.amazonaws.services.batch.model.AWSBatchException
import com.amazonaws.services.ecs.model.AccessDeniedException
import com.upplication.s3fs.S3Path
import groovy.transform.CompileDynamic
import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.cloud.aws.AmazonClientFactory
import nextflow.cloud.types.CloudMachineInfo
import nextflow.exception.AbortOperationException
import nextflow.executor.Executor
import nextflow.extension.FilesEx
import nextflow.processor.ParallelPollingMonitor
import nextflow.processor.TaskHandler
import nextflow.processor.TaskMonitor
import nextflow.processor.TaskRun
import nextflow.util.Duration
import nextflow.util.RateUnit
import nextflow.util.ServiceName
import nextflow.util.ThrottlingExecutor
/**
 * AWS Batch executor
 * https://aws.amazon.com/batch/
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@ServiceName('awsbatch')
@CompileStatic
class AwsBatchExecutor extends Executor {

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

    AwsOptions getAwsOptions() {  awsOptions  }

    /**
     * @return {@code true} to signal containers are managed directly the AWS Batch service
     */
    final boolean isContainerNative() {
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
            throw new AbortOperationException("When using `$name` executor a S3 bucket must be provided as working directory either using -bucket-dir or -work-dir command line option")
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
        def disableBinDir = session.getExecConfigProp(name, 'disableRemoteBinDir', false)
        if( session.binDir && !session.binDir.empty() && !disableBinDir ) {
            def s3 = getTempDir()
            log.info "Uploading local `bin` scripts folder to ${s3.toUriString()}/bin"
            remoteBinDir = FilesEx.copyTo(session.binDir, s3)
        }
    }

    protected void createAwsClient() {
        /*
         * retrieve config and credentials and create AWS client
         */
        final driver = new AmazonClientFactory(session.config)

        /*
         * create a proxy for the aws batch client that manages the request throttling
         */
        client = new AwsBatchProxy(driver.getBatchClient(), submitter)
        helper = createHelper(client, driver)
        // create the options object
        awsOptions = new AwsOptions(this)
        log.debug "[AWS BATCH] Executor options=$awsOptions"
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

    private AwsBatchHelper createHelper(AWSBatch batchClient, AmazonClientFactory factory) {
        try {
            return new AwsBatchHelper(
                    batchClient: batchClient,
                    ec2Client: factory.getEc2Client(),
                    ecsClient: factory.getEcsClient() )
        }
        catch (Exception e) {
            log.warn "Unable to create AWS Batch helper class | ${e.message}", e
            return null
        }
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

        final def params = [
                name: name,
                session: session,
                pollInterval: pollInterval,
                dumpInterval: dumpInterval
        ]

        log.debug "Creating parallel monitor for executor '$name' > pollInterval=$pollInterval; dumpInterval=$dumpInterval"
        new ParallelPollingMonitor(submitter, reaper, params)
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

    /**
     * @return Creates a {@link ThrottlingExecutor} service to throttle
     * the API requests to the AWS Batch service.
     */
    private ThrottlingExecutor createExecutorService(String name) {

        final qs = session.getQueueSize(name, 5_000)
        final limit = session.getExecConfigProp(name,'submitRateLimit','50/s') as String
        final size = Runtime.runtime.availableProcessors() * 5

        final opts = new ThrottlingExecutor.Options()
                .retryOn { Throwable t -> t instanceof AWSBatchException && t.errorCode=='TooManyRequestsException' }
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

    protected void logRateLimitChange(RateUnit rate) {
        log.debug "New submission rate limit: $rate"
    }

    @PackageScope
    ThrottlingExecutor getReaper() { reaper }


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

}









