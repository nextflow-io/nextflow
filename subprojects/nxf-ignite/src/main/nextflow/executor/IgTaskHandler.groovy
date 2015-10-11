/*
 * Copyright (c) 2013-2015, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2015, Paolo Di Tommaso and the respective authors.
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
import java.nio.file.NoSuchFileException
import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.exception.ProcessException
import nextflow.file.FileHelper
import nextflow.processor.TaskContext
import nextflow.processor.TaskHandler
import nextflow.processor.TaskRun
import nextflow.processor.TaskStatus
import nextflow.script.ScriptType
import nextflow.util.Duration
import org.apache.ignite.compute.ComputeJobResult
import org.apache.ignite.lang.IgniteFuture
import org.apache.ignite.lang.IgniteInClosure
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

    /**
     * The result object for this task
     */
    private IgniteFuture future

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
        final remoteTask = ( type == ScriptType.SCRIPTLET ) ? new IgScriptTask(task,sessionId) : new IgClosureTask(task,sessionId)
        future = executor.execute( remoteTask )

        future.listen( { executor.getTaskMonitor().signal(); } as IgniteInClosure )

        // mark as submitted -- transition to STARTED has to be managed by the scheduler
        status = TaskStatus.SUBMITTED
        log.trace "Task $task > Submitted"
    }

    @Override
    boolean checkIfRunning() {
        if( isSubmitted() && future ) {
            log.trace "Task ${task} > RUNNING"
            status = TaskStatus.RUNNING
            return true
        }

        return false
    }

    @Override
    boolean checkIfCompleted() {

        if( isRunning() && (future.isCancelled() || (future.isDone() && (!exitFile || readExitStatus()!=null)))  ) {
            status = TaskStatus.COMPLETED

            final result = (ComputeJobResult)future.get()
            if( result.getException() ) {
                task.error = result.getException()
                return true
            }

            if( result.isCancelled() ) {
                task.error = ProcessException.CANCELLED_ERROR
                return true
            }

            // -- the task output depend by the kind of the task executed
            if( isScriptlet() ) {
                task.stdout = outputFile
                task.stderr = errorFile
                task.exitStatus = result.getData() as Integer
            }
            else {
                def data = result.getData() as IgResultData
                task.stdout = data.value
                task.context = new TaskContext( task.processor, data.context )
            }

            log.trace "Task ${task} > DONE"
            return true
        }

        return false
    }

    protected Integer readExitStatus() {

        // read the exit status from the exit file
        final exitStatus = safeReadExitStatus()

        if( !FileHelper.workDirIsNFS ) {
            return exitStatus != null ? exitStatus : Integer.MAX_VALUE
        }

        /*
         * When working with on a NFS it may happen that the output files exists but
         * are not visible from the NFS client due to caching and network latencies
         *
         * Continue to check until a timeout is reached
         */

        if( !begin ) {
            begin = System.currentTimeMillis()
        }

        def delta = System.currentTimeMillis() - begin
        if( status == null ) {
            if( delta>timeout ) {
                // game-over
                log.warn "Unable to read command status from: $exitFile after $delta ms"
                return Integer.MAX_VALUE
            }
            log.debug "Command exit file is empty: $exitFile after $delta ms -- Try to wait a while and .. pray."
            return null
        }

        // the exit status is not null, before return check also if expected output files are available

        /*
         * When the file in a NFS folder in order to avoid false negative
         * list the content of the parent path to force refresh of NFS metadata
         * http://stackoverflow.com/questions/3833127/alternative-to-file-exists-in-java
         * http://superuser.com/questions/422061/how-to-determine-whether-a-directory-is-on-an-nfs-mounted-drive
         */
        final process = Runtime.runtime.exec("ls -la ${task.targetDir}")
        final listStatus = process.waitFor()

        if( listStatus>0 ) {
            log.debug "Can't list command target folder (exit: $listStatus): ${task.targetDir}"
            if( delta>timeout ) {
                return Integer.MAX_VALUE
            }
            return null
        }

        if( log.isTraceEnabled() ) {
            log.trace "List target folder: ${task.targetDir}\n${process.text?.indent('  ')}\n"
        }
        // destroy the process
        process.destroy()

        return exitStatus
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

        if( status ) {
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
        future?.cancel()
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