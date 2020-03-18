/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

package nextflow.executor

import java.nio.file.FileSystems
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
    TaskHandler createTaskHandler(TaskRun task) {
        assert task
        assert task.workDir

        if( task.type == ScriptType.GROOVY )
            return new NativeTaskHandler(task,this)
        else
            return new LocalTaskHandler(task,this)

    }

    @Override
    protected void register() {
        super.register()
        if( workDir.fileSystem != FileSystems.default ) {
            log.warn "Local executor only supports default file system -- Check work directory: ${getWorkDir().toUriString()}"
        }
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
    void submit() {
        // create the wrapper script
        new BashWrapperBuilder(task) .build()

        // the cmd list to launch it
        def cmd = new ArrayList(BashWrapperBuilder.BASH) << wrapperFile.getName()
        log.debug "Launch cmd line: ${cmd.join(' ')}"

        session.getExecService().submit( {

            try {
                // NOTE: make sure to redirect process output to a file otherwise
                // execution can hang when stdout/stderr is bigger than 64k
                final workDir = task.workDir.toFile()
                final logFile = new File(workDir, TaskRun.CMD_LOG)

                ProcessBuilder builder = new ProcessBuilder()
                        .redirectErrorStream(true)
                        .redirectOutput(logFile)
                        .directory(workDir)
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

        if( isSubmitted() && (process || result) ) {
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
        if( !process ) return
        log.trace("Killing process with pid: ${process.pid}")
        def cmd = "kill -TERM $process.pid"
        def proc = new ProcessBuilder('bash', '-c', cmd ).redirectErrorStream(true).start()
        def status = proc.waitFor()
        if( status != 0 ) {
            def stdout = proc.text
            log.debug "Unable to kill $task.name -- command: $cmd; exit: $status ${stdout ? "\n${stdout.indent()}":''}"
        }
    }

    /**
     * Destroy the process handler, closing all associated streams
     */
    void destroy() {

        if( destroyed ) { return }

        if( process ) {
            process.getInputStream()?.closeQuietly()
            process.getOutputStream()?.closeQuietly()
            process.getErrorStream()?.closeQuietly()
            process.destroy()
        }

        destroyed = true
    }


    /**
     * @return An {@link TraceRecord} instance holding task runtime information
     */
    @Override
    TraceRecord getTraceRecord() {
        def result = super.getTraceRecord()
        result.put('native_id', this.process?.pid)
        return result
    }

}

/**
 * Executes a native piece of groovy code
 */
@Slf4j
@CompileStatic
class NativeTaskHandler extends TaskHandler {

    Future<Object> result

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


