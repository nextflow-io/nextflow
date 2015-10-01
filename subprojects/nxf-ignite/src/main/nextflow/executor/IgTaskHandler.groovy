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

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.exception.ProcessException
import nextflow.processor.TaskContext
import nextflow.processor.TaskHandler
import nextflow.processor.TaskRun
import nextflow.processor.TaskStatus
import nextflow.script.ScriptType
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

    private IgExecutor executor

    private ScriptType type

    private Path exitFile

    private Path outputFile

    private Path errorFile

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
        final remoteTask = ( type == ScriptType.SCRIPTLET ) ? new IgBashTask(task,sessionId) : new IgClosureTask(task,sessionId)
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

        if( isRunning() && (future.isCancelled() || (future.isDone() && (!exitFile || exitFile.lastModified()>0)))  ) {
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