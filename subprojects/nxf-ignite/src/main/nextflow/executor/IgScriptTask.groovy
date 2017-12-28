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

package nextflow.executor
import java.nio.file.FileSystems
import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.file.FileHelper
import nextflow.processor.TaskRun
import org.apache.ignite.IgniteException
/**
 * Execute a remote script task into a remote Ignite cluster node
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class IgScriptTask extends IgBaseTask<Integer>   {

    private static final long serialVersionUID = -5552939711667527410L

    private transient Process process

    private transient IgScriptStagingStrategy stageStrategy

    private transient boolean isRemoteWorkDir

    private transient Path localWorkDir

    IgScriptTask( TaskRun task, UUID sessionId ) {
        super(task, sessionId)
    }

    @Override
    protected void beforeExecute() {
        isRemoteWorkDir = bean.workDir.fileSystem != FileSystems.default
        if( isRemoteWorkDir ) {
            // when work dir is allocated on a `remote` file system path
            // the input files need to be copied locally using a staging strategy
            stageStrategy = new IgScriptStagingStrategy(task: bean.clone(), sessionId: sessionId)
            stageStrategy.stage()
            // note: set staging local dir as task work dir
            localWorkDir = stageStrategy.localWorkDir
        }
        else {
            localWorkDir = bean.workDir
        }
    }

    @Override
    protected void afterExecute() {
        if( stageStrategy ) {
            stageStrategy.unstage()
            cleanupLocalWorkDir()
        }
    }

    protected void cleanupLocalWorkDir() {
        if( bean.cleanup == false ) return
        try {
            final cmd = ['bash','-c',"(sudo -n true && sudo rm -rf '$localWorkDir' || rm -rf '$localWorkDir')&>/dev/null"]
            final status = cmd.execute().waitFor()
            if( status ) log.debug "Can't cleanup path: $localWorkDir"
        }
        catch (Exception e) {
            log.debug "Error while cleaning-up path: $localWorkDir -- Cause: ${e.message ?: e}"
        }
    }

    @Override
    protected Integer execute0() throws IgniteException {

        def wrapper = new BashWrapperBuilder(bean)

        if( isRemoteWorkDir ) {
            // set custom copy strategy to turn-off wrapper script level stage/unstage input/output files
            wrapper.copyStrategy = stageStrategy
            // important: set the local work dir as script 'workDir'
            wrapper.workDir = localWorkDir
            wrapper.scratch = null
            // important add the mount for local cached files
            wrapper.containerMount = stageStrategy.localCacheDir
        }
        else {
            // set a local scratch directory if needed
            def scratch = wrapper.scratch?.toString()
            if( scratch == 'true' || scratch == 'auto' ) {
                wrapper.scratch = FileHelper.getLocalTempPath().toString()
            }
        }

        // create the shell command line to invoke
        def job = new ArrayList(bean.shell)
        job.add( wrapper.build().toString() )
        // NOTE: the actual command is wrapped by another bash whose streams
        // are redirected to null. This is important  to consume the stdout/stderr
        // of the wrapped job otherwise that output will cause the inner `tee`s hang
        List cmd = ['/bin/bash','-c', job.join(' ') + ' &>' + TaskRun.CMD_LOG]

        log.debug "Running task > ${bean.name} -- taskId=${taskId}; workdir=${localWorkDir}; remote=${isRemoteWorkDir}"
        ProcessBuilder builder = new ProcessBuilder()
                .directory(localWorkDir.toFile())
                .command(cmd)
                .redirectErrorStream(true)

        // launch and wait
        process = builder.start()
        def result = process.waitFor()

        // make sure to destroy the process and close the streams
        process.destroy()

        log.debug "Completed task > $bean.name -- taskId=${taskId}; exitStatus=$result"
        // return the exit value
        return result
    }

    @Override
    void cancel() {
        if( process ) {
            log.debug "Cancelling process for task > $bean.name -- taskId=${taskId}"
            process.destroy()
        }
        else {
            log.debug "No process to cancel for task > $bean.name -- taskId=${taskId}"
        }
    }

    String toString() {
        "${getClass().simpleName}[taskId: ${getTaskId()}; name: ${bean?.name}; workDir: $localWorkDir]"
    }
}
