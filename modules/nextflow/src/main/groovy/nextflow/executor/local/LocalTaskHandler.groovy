/*
 * Copyright 2013-2024, Seqera Labs
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
 *
 */

package nextflow.executor.local

import java.nio.file.FileSystems
import java.nio.file.Files
import java.nio.file.Path

import groovy.transform.Canonical
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.exception.ProcessException
import nextflow.exception.ProcessUnrecoverableException
import nextflow.executor.BashWrapperBuilder
import nextflow.extension.FilesEx
import nextflow.fusion.FusionAwareTask
import nextflow.fusion.FusionHelper
import nextflow.processor.TaskHandler
import nextflow.processor.TaskRun
import nextflow.processor.TaskStatus
import nextflow.trace.TraceRecord
import nextflow.util.ProcessHelper
/**
 * A process wrapper adding the ability to access to the Posix PID
 * and the {@code hasExited} flag
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class LocalTaskHandler extends TaskHandler implements FusionAwareTask {

    @Canonical
    static class TaskResult {
        Integer exitStatus
        File logs
        Throwable error
        TaskResult(int exitStatus, File logs) {
            this.logs = logs
            this.exitStatus = exitStatus
        }
        TaskResult(Throwable error) {
            this.error = error
        }
    }

    private final Path exitFile

    private final Long wallTimeMillis

    private final Path wrapperFile

    private final Path outputFile

    private final Path errorFile

    private Process process

    private boolean destroyed

    private LocalExecutor executor

    private Session session

    private volatile TaskResult result

    LocalTaskHandler(TaskRun task, LocalExecutor executor) {
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
        buildTaskWrapper()

        // create the process builder to run the task in the local computer
        final builder = createLaunchProcessBuilder()
        final logFile = builder.redirectOutput().file()
        
        // run async via thread pool
        session.getExecService().submit( {
            try {
                // start the execution and notify the event to the monitor
                process = builder.start()
                final status = process.waitFor()
                result = new TaskResult(status, logFile)
            }
            catch( Throwable ex ) {
                result = new TaskResult(ex)
            }
            finally {
                executor.getTaskMonitor().signal()
            }

        } )

        // mark as submitted -- transition to STARTED has to be managed by the scheduler
        status = TaskStatus.SUBMITTED
    }

    protected void buildTaskWrapper() {
        final wrapper = fusionEnabled()
                ? fusionLauncher()
                : new BashWrapperBuilder(task.toTaskBean())
        // create the bash command wrapper and store in the task work dir
        wrapper.build()
    }

    protected ProcessBuilder localProcessBuilder() {
        final cmd = new ArrayList<String>(BashWrapperBuilder.BASH) << wrapperFile.getName()
        log.debug "Launch cmd line: ${cmd.join(' ')}"
        // make sure it's a posix file system
        if( task.workDir.fileSystem != FileSystems.default )
            throw new ProcessUnrecoverableException("Local executor requires the use of POSIX compatible file system — offending path: ${FilesEx.toUriString(task.workDir)}")

        // NOTE: make sure to redirect process output to a file otherwise
        // execution can hang when stdout/stderr is bigger than 64k
        final workDir = task.workDir.toFile()
        final logFile = new File(workDir, TaskRun.CMD_LOG)

        return new ProcessBuilder()
                .redirectErrorStream(true)
                .redirectOutput(logFile)
                .directory(workDir)
                .command(cmd)
    }

    protected ProcessBuilder fusionProcessBuilder() {
        if( task.workDir.fileSystem == FileSystems.default )
            throw new ProcessUnrecoverableException("Fusion file system requires the use of an object storage as work directory — offending path: ${task.workDir}")

        final submit = fusionSubmitCli()
        final launcher = fusionLauncher()
        final config = task.getContainerConfig()
        final containerOpts = task.config.getContainerOptions()
        final cmd = FusionHelper.runWithContainer(launcher, config, task.getContainer(), containerOpts, submit)
        log.debug "Launch cmd line: ${cmd}"

        final logPath = Files.createTempFile('nf-task','.log')

        return new ProcessBuilder()
                .redirectErrorStream(true)
                .redirectOutput(logPath.toFile())
                .command(List.of('sh','-c', cmd))
    }

    protected ProcessBuilder createLaunchProcessBuilder() {
        return fusionEnabled()
                ? fusionProcessBuilder()
                : localProcessBuilder()
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
            task.exitStatus = result.exitStatus!=null ? result.exitStatus : Integer.MAX_VALUE
            task.error = result.error
            task.stdout = outputFile
            task.stderr = result.exitStatus && result.logs.size() ? result.logs.tail(50) : errorFile
            status = TaskStatus.COMPLETED
            destroy()
            // fusion uses a temporary file, clean it up
            if( fusionEnabled() ) result.logs.delete()
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
    protected void killTask() {
        if( !process ) return
        final pid = ProcessHelper.pid(process)
        log.trace("Killing process with pid: ${pid}")
        def cmd = "kill -TERM $pid"
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
     * @return An {@link nextflow.trace.TraceRecord} instance holding task runtime information
     */
    @Override
    TraceRecord getTraceRecord() {
        final result = super.getTraceRecord()
        if( process )
            result.put('native_id', ProcessHelper.pid(process))
        return result
    }

}
