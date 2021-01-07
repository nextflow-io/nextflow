/*
 * Copyright 2020, Microsoft Corp
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

package nextflow.cloud.azure.batch

import java.nio.file.Path

import com.azure.storage.blob.nio.AzurePath
import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import nextflow.cloud.azure.config.AzConfig
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
 * Azure Batch executor
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@ServiceName('azurebatch')
@CompileStatic
class AzBatchExecutor extends Executor implements ExtensionPoint {

    Path remoteBinDir

    AzConfig config

    AzBatchService batchService

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
        if( !(workDir instanceof AzurePath) ) {
            session.abort()
            throw new AbortOperationException("When using `$name` executor an Azure bucket must be provided as working directory either using -bucket-dir or -work-dir command line option")
        }
    }

    protected void validatePathDir() {
        def path = session.config.navigate('env.PATH')
        if( path ) {
            log.warn "Environment PATH defined in config file is ignored by Azure Batch executor"
        }
    }

    protected void uploadBinDir() {
        /*
         * upload local binaries
         */
        def disableBinDir = session.getExecConfigProp(name, 'disableRemoteBinDir', false)
        if( session.binDir && !session.binDir.empty() && !disableBinDir ) {
            final remote = getTempDir()
            log.info "Uploading local `bin` scripts folder to ${remote.toUriString()}/bin"
            remoteBinDir = FilesEx.copyTo(session.binDir, remote)
        }
    }


    /**
     * Initialise the AWS batch executor.
     */
    @Override
    protected void register() {
        super.register()
        config = AzConfig.getConfig(session)
        batchService = new AzBatchService(this)
        validateWorkDir()
        validatePathDir()
        uploadBinDir()
    }


    @PackageScope
    Path getRemoteBinDir() {
        remoteBinDir
    }

    @Override
    protected TaskMonitor createTaskMonitor() {
        TaskPollingMonitor.create(session, name, 1000, Duration.of('10 sec'))
    }

    @Override
    TaskHandler createTaskHandler(TaskRun task) {
        return new AzBatchTaskHandler(task, this)
    }

    AzBatchService getBatchService() {
        return batchService
    }

}
