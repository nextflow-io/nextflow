/*
 * Copyright 2021, Microsoft Corp
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

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.Global
import nextflow.cloud.azure.config.AzConfig
import nextflow.cloud.azure.nio.AzPath
import nextflow.exception.AbortOperationException
import nextflow.executor.Executor
import nextflow.extension.FilesEx
import nextflow.fusion.FusionHelper
import nextflow.processor.TaskHandler
import nextflow.processor.TaskMonitor
import nextflow.processor.TaskPollingMonitor
import nextflow.processor.TaskRun
import nextflow.util.Duration
import nextflow.util.ServiceName
import org.pf4j.ExtensionPoint

/**
 * Nextflow executor for Azure batch service
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@ServiceName('azurebatch')
@CompileStatic
class AzBatchExecutor extends Executor implements ExtensionPoint {

    private Path remoteBinDir

    private AzConfig config

    private AzBatchService batchService

    /**
     * @return {@code true} to signal containers are managed directly the AWS Batch service
     */
    final boolean isContainerNative() {
        return true
    }

    @Override
    String containerConfigEngine() {
        return 'docker'
    }

    @Override
    Path getWorkDir() {
        session.bucketDir ?: session.workDir
    }

    protected void validateWorkDir() {
        /*
         * make sure the work dir is an Azure bucket
         */
        if( !(workDir instanceof AzPath) ) {
            session.abort()
            throw new AbortOperationException("When using `$name` executor an Azure bucket must be provided as working directory using either the `-bucket-dir` or `-work-dir` command line option")
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
        if( session.binDir && !session.binDir.empty() && !session.disableRemoteBinDir ) {
            final remote = getTempDir()
            log.info "Uploading local `bin` scripts folder to ${remote.toUriString()}/bin"
            remoteBinDir = FilesEx.copyTo(session.binDir, remote)
        }
    }

    protected void initBatchService() {
        config = AzConfig.getConfig(session)
        batchService = new AzBatchService(this)

        // Generate an account SAS token using either activeDirectory configs or storage account keys
        if (!config.storage().sasToken) {
            config.storage().sasToken = config.activeDirectory().isConfigured() || config.managedIdentity().isConfigured()
                    ? AzHelper.generateContainerSasWithActiveDirectory(workDir, config.storage().tokenDuration)
                    : AzHelper.generateAccountSasWithAccountKey(workDir, config.storage().tokenDuration)
        }

        Global.onCleanup((it) -> batchService.close())
    }

    /**
     * Initialise the Azure Batch executor.
     */
    @Override
    protected void register() {
        super.register()
        initBatchService()
        validateWorkDir()
        validatePathDir()
        uploadBinDir()
    }

    @PackageScope AzConfig getConfig() {
        return config
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

    Path getRemoteBinDir() { return remoteBinDir }

    @Override
    boolean isFusionEnabled() {
        return FusionHelper.isFusionEnabled(session)
    }

}
