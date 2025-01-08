/*
 * Copyright 2023, Seqera Labs
 * Copyright 2022, Google Inc.
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

import com.google.cloud.storage.contrib.nio.CloudStoragePath
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.SysEnv
import nextflow.cloud.google.batch.client.BatchClient
import nextflow.cloud.google.batch.client.BatchConfig
import nextflow.cloud.google.batch.logging.BatchLogging
import nextflow.exception.AbortOperationException
import nextflow.executor.Executor
import nextflow.executor.TaskArrayExecutor
import nextflow.extension.FilesEx
import nextflow.fusion.FusionHelper
import nextflow.processor.TaskHandler
import nextflow.processor.TaskMonitor
import nextflow.processor.TaskPollingMonitor
import nextflow.processor.TaskRun
import nextflow.util.Duration
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
    private BatchConfig config
    private Path remoteBinDir
    private BatchLogging logging

    private final Set<String> deletedJobs = new HashSet<>()

    BatchClient getClient() { return client }
    BatchConfig getConfig() { return config }
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
        this.config = BatchConfig.create(session)
        log.debug "[GOOGLE BATCH] Executor config=$config"
    }

    protected void createClient() {
        this.client = new BatchClient(config)
        this.logging = new BatchLogging(config)
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
        TaskPollingMonitor.create(session, name, 1000, Duration.of('10 sec'))
    }

    @Override
    TaskHandler createTaskHandler(TaskRun task) {
        return new GoogleBatchTaskHandler(task, this)
    }

    @Override
    void shutdown() {
        client.shutdown()
        logging.close()
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
