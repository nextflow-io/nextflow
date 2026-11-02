/*
 * Copyright 2013-2026, Seqera Labs
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

package nextflow.executor.local

import java.nio.file.FileSystems
import java.nio.file.Files
import java.nio.file.Path
import java.util.concurrent.TimeUnit

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
import nextflow.util.StringUtils
import nextflow.util.SysHelper
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

    private DockerWrapperLauncher dockerLauncher

    private boolean destroyed

    private LocalExecutor executor

    private Session session

    private volatile TaskResult result

    String acceleratorEnv

    List<String> acceleratorIds

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
        executor.getExecService().submit( {
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
        if( isWindowsHost() )
            throw new ProcessUnrecoverableException(windowsErrorMessage())
        final cmd = new ArrayList<String>(BashWrapperBuilder.BASH) << wrapperFile.getName()
        log.debug "Launch cmd line: ${cmd.join(' ')}"
        // make sure it's a posix file system
        if( task.workDir.fileSystem != FileSystems.default )
            throw new ProcessUnrecoverableException("Local executor requires the use of POSIX compatible file system — offending path: ${FilesEx.toUriString(task.workDir)}")

        // NOTE: make sure to redirect process output to a file otherwise
        // execution can hang when stdout/stderr is bigger than 64k
        final workDir = task.workDir.toFile()
        final logFile = new File(workDir, TaskRun.CMD_LOG)

        final pb = new ProcessBuilder()
                .redirectErrorStream(true)
                .redirectOutput(logFile)
                .directory(workDir)
                .command(cmd)
        prepareAccelerators(pb)
        return pb
    }

    protected ProcessBuilder fusionProcessBuilder() {
        if( task.workDir.fileSystem == FileSystems.default )
            throw new ProcessUnrecoverableException("Fusion file system requires the use of an object storage as work directory — offending path: ${task.workDir}")

        final submit = fusionSubmitCli()
        final launcher = fusionLauncher()
        final config = task.getContainerConfig()
        final containerOpts = task.config.getContainerOptions()
        final cmd = FusionHelper.runWithContainer(launcher, config, task.getContainer(), containerOpts, submit, task.getContainerPlatform())
        log.debug "Launch cmd line: ${StringUtils.stripSecrets(cmd)}"

        final logPath = Files.createTempFile('nf-task','.log')

        final pb = new ProcessBuilder()
                .redirectErrorStream(true)
                .redirectOutput(logPath.toFile())
                .command(List.of('sh','-c', cmd))
        prepareAccelerators(pb)
        return pb
    }

    /**
     * Run the task wrapper script inside the Docker container, see {@link LocalExecutor#isDockerWrapperEnabled()}
     */
    protected ProcessBuilder dockerProcessBuilder() {
        // make sure it's a local file system, so it can be mounted in the container
        if( task.workDir.fileSystem != FileSystems.default )
            throw new ProcessUnrecoverableException("Local executor requires the use of a local file system — offending path: ${FilesEx.toUriString(task.workDir)}")

        dockerLauncher = new DockerWrapperLauncher(task.toTaskBean(), containerName())
        final cmd = dockerLauncher.runCommand()
        log.debug "Launch cmd line: ${StringUtils.stripSecrets(cmd.join(' '))}"

        // NOTE: make sure to redirect process output to a file otherwise
        // execution can hang when stdout/stderr is bigger than 64k
        final workDir = task.workDir.toFile()
        final logFile = new File(workDir, TaskRun.CMD_LOG)

        final pb = new ProcessBuilder()
                .redirectErrorStream(true)
                .redirectOutput(logFile)
                .directory(workDir)
                .command(cmd)
        prepareAccelerators(pb)
        return pb
    }

    protected String containerName() {
        final hash = task.hash?.toString()
        final prefix = hash ? hash.substring(0, Math.min(12, hash.length())) + '-' : ''
        return 'nxf-' + prefix + UUID.randomUUID().toString().substring(0, 8)
    }

    protected boolean isWindowsHost() {
        return SysHelper.isWindows()
    }

    protected String windowsErrorMessage() {
        if( executor.isDockerWrapperEnabled() )
            return "Process `${task.lazyName()}` does not specify a container -- On Windows, every process must define a container to be executed with Docker"
        return "Process `${task.lazyName()}` cannot be executed -- On Windows, tasks can only be executed with Docker; enable it with the setting `docker.enabled = true` and define a container for each process"
    }

    protected void prepareAccelerators(ProcessBuilder pb) {
        if( !acceleratorEnv )
            return
        pb.environment().put(acceleratorEnv, acceleratorIds.join(','))
    }

    protected ProcessBuilder createLaunchProcessBuilder() {
        if( fusionEnabled() )
            return fusionProcessBuilder()
        if( executor.isDockerWrapperEnabled() && task.getContainer() )
            return dockerProcessBuilder()
        return localProcessBuilder()
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
                // destroy() is async, so wait for the process to actually exit
                task.exitStatus = process.waitFor(5, TimeUnit.SECONDS) ? process.exitValue() : Integer.MAX_VALUE
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
        if( dockerLauncher ) {
            stopContainer()
            return
        }
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
     * Stop the task container. Killing the `docker run` client process is not enough
     * because the container keeps running (on Windows, signals are not forwarded)
     */
    protected void stopContainer() {
        final cmd = dockerLauncher.stopCommand()
        if( !cmd ) {
            process.destroy()
            return
        }
        log.trace "Stopping task container: ${cmd.join(' ')}"
        try {
            // note: do not wait for the command completion, the container is stopped asynchronously
            new ProcessBuilder(cmd)
                .redirectErrorStream(true)
                .redirectOutput(ProcessBuilder.Redirect.DISCARD)
                .start()
        }
        catch( Exception e ) {
            log.debug "Unable to stop container for task $task.name -- command: ${cmd.join(' ')}", e
        }
    }

    /**
     * Destroy the process handler, closing all associated streams
     */
    void destroy() {

        if( destroyed ) { return }

        // stop the container when the task is still running e.g. it exceeded the time limit
        if( dockerLauncher && process?.isAlive() )
            stopContainer()

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
        if( dockerLauncher )
            result.put('native_id', dockerLauncher.getName())
        else if( process )
            result.put('native_id', ProcessHelper.pid(process))
        return result
    }

}
