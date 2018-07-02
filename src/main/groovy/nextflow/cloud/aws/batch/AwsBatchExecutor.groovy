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

import com.amazonaws.services.batch.AWSBatchClient
import com.upplication.s3fs.S3Path
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
     * AWS batch client instance
     */
    @PackageScope
    private static AWSBatchClient client

    private static Path remoteBinDir = null

    final boolean isContainerNative() {
        return true
    }

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
        client = new AmazonCloudDriver(session.config).getBatchClient()

    }

    @PackageScope
    Path getRemoteBinDir() {
        remoteBinDir
    }

    @PackageScope
    AWSBatchClient getClient() {
        client
    }

    @Override
    protected TaskMonitor createTaskMonitor() {
        ParallelPollingMonitor.create(session, name, Duration.of('10 sec'))
    }

    @Override
    TaskHandler createTaskHandler(TaskRun task) {
        assert task
        assert task.workDir
        log.trace "[AWS BATCH] Launching process > ${task.name} -- work folder: ${task.workDirStr}"
        new AwsBatchTaskHandler(task, this)
    }
}









