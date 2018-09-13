/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
 */

package nextflow.cloud.aws.batch

import java.nio.file.Path

import com.amazonaws.services.batch.AWSBatch
import com.amazonaws.services.batch.model.AWSBatchException
import com.upplication.s3fs.S3Path
import groovy.transform.CompileDynamic
import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.Nextflow
import nextflow.cloud.aws.AmazonCloudDriver
import nextflow.exception.AbortOperationException
import nextflow.executor.Executor
import nextflow.extension.FilesEx
import nextflow.processor.ParallelPollingMonitor
import nextflow.processor.TaskHandler
import nextflow.processor.TaskMonitor
import nextflow.processor.TaskRun
import nextflow.util.Duration
import nextflow.util.RateUnit
import nextflow.util.ThrottlingExecutor
/**
 * AWS Batch executor
 * https://aws.amazon.com/batch/
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class AwsBatchExecutor extends Executor {

    /**
     * Proxy to throttle AWS batch client requests
     */
    @PackageScope
    private static AwsBatchProxy client

    /**
     * executor service to throttle service requests
     */
    private static ThrottlingExecutor submitter

    /**
     * Executor service to throttle cancel requests
     */
    private static ThrottlingExecutor reaper

    /**
     * A S3 path where executable scripts need to be uploaded
     */
    private static Path remoteBinDir = null

    /**
     * @return {@code true} to signal containers are managed directly the AWS Batch service
     */
    final boolean isContainerNative() {
        return true
    }

    /**
     * Initialise the AWS batch executor.
     */
    @Override
    void register() {
        super.register()

        /*
         * make sure the work dir is a S3 bucket
         */
        if( !(session.workDir instanceof S3Path) ) {
            session.abort()
            throw new AbortOperationException("When using `$name` executor a S3 bucket must be provided as working directory -- Add the option `-w s3://<your-bucket/path>` to your run command line")
        }

        def path = session.config.navigate('env.PATH')
        if( path ) {
            log.warn "Environment PATH defined in config file is ignored by AWS Batch executor"
        }

        /*
         * upload local binaries
         */
        def disableBinDir = session.getExecConfigProp(name, 'disableRemoteBinDir', false)
        if( session.binDir && !session.binDir.empty() && !disableBinDir ) {
            def s3 = Nextflow.tempDir()
            log.info "Uploading local `bin` scripts folder to ${s3.toUriString()}/bin"
            remoteBinDir = FilesEx.copyTo(session.binDir, s3)
        }

        /*
         * retrieve config and credentials and create AWS client
         */
        final driver = new AmazonCloudDriver(session.config)

        /*
         * create a proxy for the aws batch client that manages the request throttling 
         */
        client = new AwsBatchProxy(driver.getBatchClient(), submitter)
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

}









