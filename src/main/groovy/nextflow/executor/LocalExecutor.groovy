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
import java.nio.file.Path
import java.util.concurrent.Callable
import java.util.concurrent.Future

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.exception.ProcessException
import nextflow.processor.LocalPollingMonitor
import nextflow.processor.TaskHandler
import nextflow.processor.TaskMonitor
import nextflow.processor.TaskRun
import nextflow.processor.TaskStatus
import nextflow.script.ScriptType
import nextflow.trace.TraceRecord
import nextflow.util.Escape
import nextflow.util.PosixProcess
/**
 * Executes the specified task on the locally exploiting the underlying Java thread pool
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
@SupportedScriptTypes( [ScriptType.SCRIPTLET, ScriptType.GROOVY] )
class LocalExecutor extends Executor {

    @Override
    protected TaskMonitor createTaskMonitor() {

        return LocalPollingMonitor.create(session, name)

    }

    @Override
    def TaskHandler createTaskHandler(TaskRun task) {
        assert task
        assert task.workDir

        if( task.type == ScriptType.GROOVY )
            return new NativeTaskHandler(task,this)
        else
            return new LocalTaskHandler(task,this)

    }

}


/**
 * A process wrapper adding the ability to access to the Posix PID
 * and the {@code hasExited} flag
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class LocalTaskHandler extends TaskHandler {

    private final Path exitFile

    private final Long wallTimeMillis

    private final Path wrapperFile

    private final Path outputFile

    private final Path errorFile

    private PosixProcess process

    private boolean destroyed

    private LocalExecutor executor

    private Session session

    private volatile result


    LocalTaskHandler( TaskRun task, LocalExecutor executor  ) {
        super(task)
        // create the task handler
        this.exitFile = task.workDir.resolve(TaskRun.CMD_EXIT)
        this.outputFile = task.workDir.resolve(TaskRun.CMD_OUTFILE)
        this.errorFile = task.workDir.resolve(TaskRun.CMD_ERRFILE)
        this.wrapperFile = task.workDir.resolve(TaskRun.CMD_RUN)
        this.wallTimeMillis = task.config.getTime()?.toMillis()
        this.executor = executor
        this.session = executor.session
    }


    @Override
    def void submit() {
        // create the wrapper script
        new BashWrapperBuilder(task) .build()

        // the cmd list to launch it
        def job = new ArrayList(BashWrapperBuilder.BASH) << wrapperFile.getName()
        // NOTE: the actual command is wrapper by another bash whose stream
        // are redirected to null. This is important in order to consume the stdout/stderr
        // of the wrapped job otherwise that output will cause the inner `tee`s hang
        List cmd = ['/bin/bash','-c', job.join(' ') + " &> ${Escape.path(task.workDir)}/${TaskRun.CMD_LOG}" ]
        log.trace "Launch cmd line: ${cmd.join(' ')}"

        session.getExecService().submit( {

            try {
                ProcessBuilder builder = new ProcessBuilder()
                        .directory(task.workDir.toFile())
                        .command(cmd)

                // -- start the execution and notify the event to the monitor
                process = new PosixProcess(builder.start())
                result = process.waitFor()
            }
            catch( Throwable ex ) {
                log.trace("Failed to execute command: ${cmd.join(' ')}", ex)
                result = ex
            }
            finally {
                executor.getTaskMonitor().signal()
            }

        } )

        // mark as submitted -- transition to STARTED has to be managed by the scheduler
        status = TaskStatus.SUBMITTED
    }


    long elapsedTimeMillis() {
        startTimeMillis ? System.currentTimeMillis() - startTimeMillis : 0
    }

    /**
     * Check if the submitted job has started
     */
    @Override
    boolean checkIfRunning() {

        if( isSubmitted() && process != null ) {
            status = TaskStatus.RUNNING
            return true
        }

        return false
    }

    /**
     * Check if the submitted job has terminated its execution
     */
    @Override
    boolean checkIfCompleted() {

        if( !isRunning() ) { return false }

        if( result != null ) {
            task.exitStatus = result instanceof Integer ? result : Integer.MAX_VALUE
            task.error = result instanceof Throwable ? result : null
            task.stdout = outputFile
            task.stderr = errorFile
            status = TaskStatus.COMPLETED
            destroy()
            return true
        }

        if( wallTimeMillis ) {
            /*
             * check if the task exceed max duration time
             */
            if( elapsedTimeMillis() > wallTimeMillis ) {
                destroy()
                task.stdout = outputFile
                task.stderr = errorFile
                task.error = new ProcessException("Process exceeded running time limit (${task.config.getTime()})")
                status = TaskStatus.COMPLETED

                // signal it has completed
                return true
            }
        }

        return false
    }


    /**
     * Force the submitted job to quit
     */
    @Override
    void kill() {
        log.trace("Killing process with pid: ${process.pid}")
        def pid = process.pid.toString()
        new ProcessBuilder('bash', '-c', "kill -TERM $pid" ).start().waitFor()
    }

    /**
     * Destroy the process handler, closing all associated streams
     */
    void destroy() {

        if( destroyed ) { return }

        process.getInputStream()?.closeQuietly()
        process.getOutputStream()?.closeQuietly()
        process.getErrorStream()?.closeQuietly()
        process.destroy()
        destroyed = true
    }


    /**
     * @return An {@link TraceRecord} instance holding task runtime information
     */
    @Override
    TraceRecord getTraceRecord() {
        def result = super.getTraceRecord()
        result.native_id = this.process?.pid
        return result
    }

}

/**
 * Executes a native piece of groovy code
 */
@Slf4j
@CompileStatic
class NativeTaskHandler extends TaskHandler {

    def Future<Object> result

    private Session session

    private Executor executor

    private class TaskSubmit implements Callable {

        final TaskRun task

        TaskSubmit( TaskRun obj ) { task = obj }

        @Override
        Object call() throws Exception {
            try  {
                return task.code.call()
            }
            catch( Throwable error ) {
                return error
            }
            finally {
                executor.getTaskMonitor().signal()
            }
        }
    }

    protected NativeTaskHandler(TaskRun task, Executor executor) {
        super(task)
        this.executor = executor
        this.session = executor.session
    }


    @Override
    void submit() {
        // submit for execution by using session executor service
        // it returns an error when everything is OK
        // of the exception throw in case of error
        result = session.getExecService().submit(new TaskSubmit(task))
        status = TaskStatus.SUBMITTED
    }

    @Override
    boolean checkIfRunning() {
        if( isSubmitted() && result != null ) {
            status = TaskStatus.RUNNING
            return true
        }

        return false
    }

    @Override
    boolean checkIfCompleted() {
        if( isRunning() && result.isDone() ) {
            status = TaskStatus.COMPLETED
            if( result.get() instanceof Throwable ) {
                task.error = (Throwable)result.get()
            }
            else {
                task.stdout = result.get()
            }
            return true
        }
        return false
    }

    @Override
    void kill() {
        if( result ) result.cancel(true)
    }

}


