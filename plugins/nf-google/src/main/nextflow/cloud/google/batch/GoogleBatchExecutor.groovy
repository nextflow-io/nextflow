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

package nextflow.cloud.google.batch

import static nextflow.cloud.google.batch.GoogleBatchScriptLauncher.*

import java.nio.file.Path
import java.util.concurrent.TimeUnit

import com.google.api.gax.rpc.DeadlineExceededException
import com.google.api.gax.rpc.UnavailableException
import com.google.cloud.storage.contrib.nio.CloudStoragePath
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.SysEnv
import nextflow.cloud.google.GoogleOpts
import nextflow.cloud.google.batch.client.BatchClient
import nextflow.cloud.google.batch.client.BatchConfig
import nextflow.cloud.google.batch.logging.BatchLogging
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
import nextflow.util.RateUnit
import nextflow.util.ThreadPoolHelper
import nextflow.util.ThrottlingExecutor
import nextflow.util.Escape
import nextflow.util.ServiceName
import org.pf4j.ExtensionPoint
/**
 * Implements support for Google Batch
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@ServiceName(value='google-batch')
@CompileStatic
class GoogleBatchExecutor extends Executor implements ExtensionPoint, TaskArrayExecutor {

    private BatchClient client
    private GoogleOpts googleOpts
    private Path remoteBinDir
    private BatchLogging logging

    /**
     * Executor service to throttle job submission requests
     */
    private ThrottlingExecutor submitter

    /**
     * Executor service to throttle job deletion requests
     */
    private ThrottlingExecutor reaper

    private final Set<String> deletedJobs = new HashSet<>()

    BatchClient getClient() { return client }
    GoogleOpts getGoogleOpts() { return googleOpts }
    BatchConfig getBatchConfig() { return googleOpts.batch }
    Path getRemoteBinDir() { return remoteBinDir }
    BatchLogging getLogging() { logging }

    @Override
    final boolean isContainerNative() {
        return true
    }

    @Override
    String containerConfigEngine() {
        return 'docker'
    }

    @Override
    final Path getWorkDir() {
        return session.bucketDir ?: session.workDir
    }

    protected void validateWorkDir() {
        if ( getWorkDir()?.scheme != 'gs' ) {
            session.abort()
            throw new AbortOperationException("Executor `google-batch` requires a Google Storage bucket to be specified as a working directory -- Add the option `-w gs://<your-bucket/path>` to your run command line or specify a workDir in your config file")
        }
    }

    protected void uploadBinDir() {
        if( session.binDir && !session.binDir.empty() && !session.disableRemoteBinDir ) {
            final cloudPath = getTempDir()
            log.info "Uploading local `bin` scripts folder to ${cloudPath.toUriString()}/bin"
            this.remoteBinDir = FilesEx.copyTo(session.binDir, cloudPath)
        }
    }

    protected void createConfig() {
        this.googleOpts = GoogleOpts.create(session)
        log.debug "[GOOGLE BATCH] Executor config=$googleOpts"
    }

    protected void createClient() {
        this.client = new BatchClient(googleOpts)
        this.logging = new BatchLogging(googleOpts)
    }

    @Override
    protected void register() {
        super.register()
        createConfig()
        validateWorkDir()
        uploadBinDir()
        createClient()
    }

    @Override
    protected TaskMonitor createTaskMonitor() {
        // create the throttling executor services
        submitter = createExecutorService('GoogleBatch-executor')
        reaper = createExecutorService('GoogleBatch-reaper')

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

        log.debug "Creating parallel monitor for executor '$name' > capacity=$capacity; pollInterval=$pollInterval; dumpInterval=$dumpInterval"
        new ParallelPollingMonitor(submitter, params)
    }

    /**
     * Creates a {@link ThrottlingExecutor} service to throttle
     * API requests to the Google Batch service.
     *
     * @param name The executor service name
     * @return A {@link ThrottlingExecutor} instance
     */
    private ThrottlingExecutor createExecutorService(String name) {
        final qs = 5_000
        final limit = config.getExecConfigProp(name, 'submitRateLimit', '50/s') as String
        final size = Runtime.runtime.availableProcessors() * 5

        final opts = new ThrottlingExecutor.Options()
                .retryOn { Throwable t ->
                    t instanceof UnavailableException ||
                    t instanceof DeadlineExceededException ||
                    t instanceof IOException ||
                    t.cause instanceof IOException
                }
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

    protected void logRateLimitChange(RateUnit rate) {
        log.debug "New submission rate limit: $rate"
    }

    @Override
    TaskHandler createTaskHandler(TaskRun task) {
        return new GoogleBatchTaskHandler(task, this)
    }

    ThrottlingExecutor getReaper() { reaper }

    @Override
    void shutdown() {
        // shutdown the submitter executor
        def tasks = submitter.shutdownNow()
        if( tasks ) log.warn "Execution interrupted -- cleaning up execution pool"
        submitter.awaitTermination(5, TimeUnit.MINUTES)

        // shutdown the reaper executor
        reaper.shutdown()
        final waitMsg = "[GOOGLE BATCH] Waiting jobs reaper to complete (%d jobs to be terminated)"
        final exitMsg = "[GOOGLE BATCH] Exiting before jobs reaper thread pool complete -- Some jobs may not be terminated"
        awaitCompletion(reaper, Duration.of('60min'), waitMsg, exitMsg)

        // close batch client and logging
        client.shutdown()
        logging.close()
    }

    protected void awaitCompletion(ThrottlingExecutor executor, Duration duration, String waitMsg, String exitMsg) {
        try {
            ThreadPoolHelper.await(executor, duration, waitMsg, exitMsg)
        }
        catch( java.util.concurrent.TimeoutException e ) {
            log.warn(e.message, e)
        }
    }

    @Override
    boolean isFusionEnabled() {
        return FusionHelper.isFusionEnabled(session)
    }

    boolean isCloudinfoEnabled() {
        return Boolean.parseBoolean(SysEnv.get('NXF_CLOUDINFO_ENABLED', 'true') )
    }

    boolean shouldDeleteJob(String jobId) {
        if( jobId in deletedJobs ) {
            // if the job is already in the list it has been already deleted
            return false
        }
        synchronized (deletedJobs) {
            // add the job id to the set of deleted jobs, if it's a new id, the `add` method
            // returns true therefore the job should be deleted
            return deletedJobs.add(jobId)
        }
    }

    @Override
    String getArrayIndexName() {
        return 'BATCH_TASK_INDEX'
    }

    @Override
    int getArrayIndexStart() {
        return 0
    }

    @Override
    String getArrayTaskId(String jobId, int index) {
        return index.toString()
    }

    @Override
    String getArrayWorkDir(TaskHandler handler) {
        return isFusionEnabled() || isWorkDirDefaultFS()
            ? TaskArrayExecutor.super.getArrayWorkDir(handler)
            : containerMountPath(handler.task.workDir as CloudStoragePath)
    }

    @Override
    String getArrayLaunchCommand(String taskDir) {
        if( isFusionEnabled() || isWorkDirDefaultFS() ) {
            return TaskArrayExecutor.super.getArrayLaunchCommand(taskDir)
        }
        else {
            final cmd = List.of('/bin/bash','-o','pipefail','-c', launchCommand(taskDir))
            return Escape.cli(cmd as String[])
        }
    }

}
