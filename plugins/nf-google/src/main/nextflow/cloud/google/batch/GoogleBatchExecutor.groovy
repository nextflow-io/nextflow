/*
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

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.cloud.google.batch.client.BatchConfig
import nextflow.cloud.google.batch.client.BatchClient
import nextflow.cloud.google.batch.logging.BatchLogging
import nextflow.exception.AbortOperationException
import nextflow.executor.Executor
import nextflow.extension.FilesEx
import nextflow.processor.TaskHandler
import nextflow.processor.TaskMonitor
import nextflow.processor.TaskPollingMonitor
import nextflow.processor.TaskRun
import nextflow.util.Duration
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
class GoogleBatchExecutor extends Executor implements ExtensionPoint {

    private BatchClient client
    private BatchConfig config
    private Path remoteBinDir
    private BatchLogging logging

    BatchClient getClient() { return client }
    BatchConfig getConfig() { return config }
    Path getRemoteBinDir() { return remoteBinDir }
    BatchLogging getLogging() { logging }

    @Override
    final boolean isContainerNative() {
        return true
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
    }
}
