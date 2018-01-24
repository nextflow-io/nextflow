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
import static nextflow.processor.TaskStatus.COMPLETED
import static nextflow.processor.TaskStatus.RUNNING
import static nextflow.processor.TaskStatus.SUBMITTED

import java.nio.file.NoSuchFileException
import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.file.FileHelper
import nextflow.processor.TaskContext
import nextflow.processor.TaskHandler
import nextflow.processor.TaskRun
import nextflow.script.ScriptType
import nextflow.util.Duration
/**
 * A task handler for Ignite  cluster
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class IgTaskHandler extends TaskHandler {

    final static private Duration READ_TIMEOUT = Duration.of('270sec') // 4.5 minutes

    private long timeout

    private IgExecutor executor

    private ScriptType type

    private Path exitFile

    private Path outputFile

    private Path errorFile

    private long begin

    static IgTaskHandler createScriptHandler( TaskRun task, IgExecutor executor ) {
        def handler = new IgTaskHandler(task)
        handler.executor = executor
        handler.type = ScriptType.SCRIPTLET
        handler.exitFile = task.workDir.resolve(TaskRun.CMD_EXIT)
        handler.outputFile = task.workDir.resolve(TaskRun.CMD_OUTFILE)
        handler.errorFile = task.workDir.resolve(TaskRun.CMD_ERRFILE)
        handler.timeout = (executor.session?.getExitReadTimeout(executor.name, READ_TIMEOUT) ?: READ_TIMEOUT).toMillis()
        return handler
    }

    static IgTaskHandler createGroovyHandler( TaskRun task, IgExecutor executor ) {
        def handler = new IgTaskHandler(task)
        handler.executor = executor
        handler.type = ScriptType.GROOVY
        return handler
    }

    private IgTaskHandler(TaskRun task) {
        super(task)
    }

    @Override
    void submit() {

        // submit to a Ignite node for execution
        final sessionId = task.processor.session.uniqueId
        final remoteTask = type == ScriptType.SCRIPTLET ? new IgScriptTask(task,sessionId) : new IgClosureTask(task,sessionId)
        executor.execute( remoteTask )

        // mark as submitted -- transition to STARTED has to be managed by the scheduler
        status = SUBMITTED
        log.trace "Task SUBMITTED > $task"
    }

    @Override
    boolean checkIfRunning() {
        if( isSubmitted() && isStarted() ) {
            log.trace "Task RUNNING > ${task}"
            status = RUNNING
            return true
        }

        return false
    }

    private boolean isStarted() {
        executor.checkTaskStarted(task.id)
    }

    @Override
    boolean checkIfCompleted() {

        if( isRunning() && executor.checkTaskCompleted(task.id) ) {

            if( !executor.checkTaskFailed(task.id) && exitFile && readExitStatus()==null ) {
                // continue to wait for a valid exit status
                return false
            }

            log.trace "Task DONE > ${task}"
            status = COMPLETED
            final holder = executor.removeTaskCompleted(task.id)

            if( holder.error ) {
                task.error = holder.error
                return true
            }

            // -- the task output depend by the kind of the task executed
            if( isScriptlet() ) {
                task.stdout = outputFile
                task.stderr = errorFile
                task.exitStatus = holder.result as Integer
            }
            else {
                def data = holder.result as IgResultData
                task.stdout = data.value
                task.context = new TaskContext( task.processor, data.context )
            }
            return true
        }

        return false
    }

    protected Integer readExitStatus() {

        String workDirList = null
        if( !begin ) {
            begin = System.currentTimeMillis()
        }
        else {
            /*
             * When the file is in a NFS folder in order to avoid false negative
             * list the content of the parent path to force refresh of NFS metadata
             * http://stackoverflow.com/questions/3833127/alternative-to-file-exists-in-java
             * http://superuser.com/questions/422061/how-to-determine-whether-a-directory-is-on-an-nfs-mounted-drive
             */
            workDirList = FileHelper.listDirectory(task.workDir)
        }

        // read the exit status from the exit file
        final exitStatus = safeReadExitStatus()
        if( exitStatus != null )
            return exitStatus

        if( !FileHelper.workDirIsNFS ) {
            return Integer.MAX_VALUE
        }

        /*
         * When working with on a NFS it may happen that the output files exists but
         * are not visible from the NFS client due to caching and network latencies
         *
         * Continue to check until a timeout is reached
         */

        def delta = System.currentTimeMillis() - begin
        if( delta > timeout ) {
            // game-over
            def errMessage = []
            errMessage << "Failed to get exist status for process ${this} -- exit-status-read-timeout=${Duration.of(timeout)}; delta=${Duration.of(delta)}"
            // -- dump current queue stats
            errMessage << "Current queue status:"
            errMessage << executor.dumpQueueStatus()?.indent('> ')
            // -- dump directory listing
            errMessage << "Content of workDir: ${task.workDir}"
            errMessage << workDirList?.indent('> ')
            log.debug errMessage.join('\n')

            return Integer.MAX_VALUE
        }

        log.debug "Command exit file is empty: $exitFile after ${Duration.of(delta)} -- Try to wait a while and .. pray."
        return null
    }


    private Integer safeReadExitStatus() {

        /*
         * read the exit file, it should contain the executed process exit status
         */
        def status

        try {
            status = exitFile.text
        }
        catch( NoSuchFileException e )  {
            return null
        }

        if( status != null ) {
            try {
                return status.toInteger()
            }
            catch( NumberFormatException e ) {
                log.warn "Unable to parse process exit file: $exitFile -- bad value: '$status'"
            }
        }

        return null
    }


    @Override
    void kill() {
        executor.cancelTask(task.id)
    }

    /**
     * @return Whenever is a shell script task
     */
    boolean isScriptlet() { type == ScriptType.SCRIPTLET }

    /**
     * @return Whenever is a groovy closure task
     */
    boolean isGroovy() { type == ScriptType.GROOVY }

}