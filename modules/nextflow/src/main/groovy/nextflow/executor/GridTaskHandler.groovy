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

package nextflow.executor

import static nextflow.processor.TaskStatus.*

import java.nio.file.Path
import java.nio.file.attribute.BasicFileAttributes
import java.time.temporal.ChronoUnit
import java.util.regex.Pattern

import dev.failsafe.Failsafe
import dev.failsafe.RetryPolicy
import dev.failsafe.event.EventListener
import dev.failsafe.event.ExecutionAttemptedEvent
import dev.failsafe.function.CheckedPredicate
import dev.failsafe.function.CheckedSupplier
import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.util.logging.Slf4j
import nextflow.Global
import nextflow.exception.ProcessException
import nextflow.exception.ProcessFailedException
import nextflow.exception.ProcessNonZeroExitStatusException
import nextflow.file.FileHelper
import nextflow.fusion.FusionAwareTask
import nextflow.fusion.FusionHelper
import nextflow.processor.TaskArrayRun
import nextflow.processor.TaskHandler
import nextflow.processor.TaskRun
import nextflow.trace.TraceRecord
import nextflow.util.CmdLineHelper
import nextflow.util.Duration
import nextflow.util.TestOnly
import nextflow.util.Throttle
/**
 * Handles a job execution in the underlying grid platform
 */
@Slf4j
@CompileStatic
class GridTaskHandler extends TaskHandler implements FusionAwareTask {

    /** The target executor platform */
    final AbstractGridExecutor executor

    /** Location of the file created when the job is started */
    final Path startFile

    /** Location of the file created when the job is terminated */
    final Path exitFile

    /** Location of the file holding the task std output */
    final Path outputFile

    /** Location of the file holding the task std error */
    final Path errorFile

    /** The wrapper file used to execute the user script */
    final Path wrapperFile

    /** The unique job ID as provided by the underlying grid platform */
    private jobId

    private queue

    private long exitStatusReadTimeoutMillis

    private Duration sanityCheckInterval

    BatchCleanup batch

    @TestOnly
    protected GridTaskHandler() {
        this.exitAwaiter = new ExitStatusAwaiter(Duration.of('270sec'))
    }

    GridTaskHandler( TaskRun task, AbstractGridExecutor executor ) {
        super(task)

        this.executor = executor
        this.startFile = task.workDir.resolve(TaskRun.CMD_START)
        this.exitFile = task.workDir.resolve(TaskRun.CMD_EXIT)
        this.outputFile = task.workDir.resolve(TaskRun.CMD_OUTFILE)
        this.errorFile = task.workDir.resolve(TaskRun.CMD_ERRFILE)
        this.wrapperFile = task.workDir.resolve(TaskRun.CMD_RUN)
        final duration = executor.config.getExitReadTimeout(executor.name)
        this.exitStatusReadTimeoutMillis = duration.toMillis()
        this.exitAwaiter = new ExitStatusAwaiter(duration)
        this.queue = task.config?.queue
        this.sanityCheckInterval = duration
    }

    @Override
    void prepareLauncher() {
        // -- create the wrapper script
        createTaskWrapper(task).build()
    }

    protected ProcessBuilder createProcessBuilder() {

        // -- log the qsub command
        final cli = executor.getSubmitCommandLine(task, wrapperFile)
        log.trace "start process ${task.name} > cli: ${cli}"

        /*
         * launch 'sub' script wrapper
         */
        ProcessBuilder builder = new ProcessBuilder()
            .command( cli as String[] )
            .redirectErrorStream(true)
        if( !fusionEnabled() )
            builder .directory(task.workDir.toFile())

        return builder
    }

    @Memoized
    protected CheckedPredicate<? extends Throwable> retryCondition(String reasonPattern) {
        final pattern = Pattern.compile(reasonPattern)
        return new CheckedPredicate<Throwable>() {
            @Override
            boolean test(Throwable failure) {
                if( failure instanceof ProcessNonZeroExitStatusException ) {
                    final reason = failure.reason
                    return reason ? pattern.matcher(reason).find() : false
                }
                return false
            }
        }
    }

    protected <T> RetryPolicy<T> retryPolicy() {

        final retry = executor.config.retry
        final listener = new EventListener<ExecutionAttemptedEvent>() {
            @Override
            void accept(ExecutionAttemptedEvent event) throws Throwable {
                final failure = event.getLastException()
                if( failure instanceof ProcessNonZeroExitStatusException ) {
                    final failure0 = (ProcessNonZeroExitStatusException)failure
                    final msg = """\
                        Failed to submit process '${task.name}'
                         - attempt : ${event.attemptCount}
                         - command : ${failure0.command}
                         - reason  : ${failure0.reason}
                        """.stripIndent(true)
                    log.warn msg

                } else {
                    log.debug("Unexpected retry failure: ${failure?.message}", failure)
                }
            }
        }

        return RetryPolicy.<T>builder()
                .handleIf(retryCondition(retry.reason))
                .withBackoff(retry.delay.toMillis(), retry.maxDelay.toMillis(), ChronoUnit.MILLIS)
                .withMaxAttempts(retry.maxAttempts)
                .withJitter(retry.jitter)
                .onFailedAttempt(listener)
                .build()
    }

    protected <T> T safeExecute(CheckedSupplier<T> action) {
        final policy = retryPolicy()
        return Failsafe.with(policy).get(action)
    }

    protected String processStart(ProcessBuilder builder, String pipeScript) {
        final process = builder.start()

        try {
            // -- forward the job launcher script to the command stdin if required
            if( pipeScript ) {
                log.trace "[${executor.name.toUpperCase()}] Submit STDIN command ${task.name} >\n${pipeScript.indent()}"
                process.out << pipeScript
                process.out.close()
            }

            // -- wait the the process completes
            final result = process.text
            final exitStatus = process.waitFor()
            final cmd = launchCmd0(builder,pipeScript)

            if( exitStatus ) {
                throw new ProcessNonZeroExitStatusException("Failed to submit process to grid scheduler for execution", result, exitStatus, cmd)
            }

            // -- return the process stdout
            return result
        }
        finally {
            // make sure to release all resources
            process.in.closeQuietly()
            process.out.closeQuietly()
            process.err.closeQuietly()
            process.destroy()
        }
    }

    protected BashWrapperBuilder createTaskWrapper(TaskRun task) {
        return fusionEnabled()
            ? fusionLauncher()
            : executor.createBashWrapperBuilder(task)
    }

    protected String stdinLauncherScript() {
        return fusionEnabled() ? fusionStdinWrapper() : wrapperFile.text
    }

    protected String fusionStdinWrapper() {
        final submit = fusionSubmitCli()
        final launcher = fusionLauncher()
        final config = task.getContainerConfig()
        final containerOpts = task.config.getContainerOptions()
        final cmd = FusionHelper.runWithContainer(launcher, config, task.getContainer(), containerOpts, submit, task.getContainerPlatform())
        // create an inline script to launch the job execution
        return '#!/bin/bash\n' + submitDirective(task) + cmd + '\n'
    }

    protected String submitDirective(TaskRun task) {
        final remoteLog = task.workDir.resolve(TaskRun.CMD_LOG).toString()
        // replaces the log file with a null file because the cluster submit tool
        // cannot write to a file hosted in a remote object storage
        final result = executor
                .getHeaders(task)
                .replaceAll(remoteLog, '/dev/null')
        return result
    }

    protected String launchCmd0(ProcessBuilder builder, String pipeScript) {
        def result = CmdLineHelper.toLine(builder.command())
        if( pipeScript ) {
            result = "cat << 'LAUNCH_COMMAND_EOF' | ${result}\n"
            result += pipeScript.trim() + '\n'
            result += 'LAUNCH_COMMAND_EOF\n'
        }
        return result
    }

    /*
     * {@inheritDocs}
     */
    @Override
    void submit() {
        ProcessBuilder builder = null
        try {
            // -- start the execution and notify the event to the monitor
            builder = createProcessBuilder()
            // -- forward the job launcher script to the command stdin if required
            final stdinScript = executor.pipeLauncherScript() ? stdinLauncherScript() : null
            // -- execute with a re-triable strategy
            final result = safeExecute( () -> processStart(builder, stdinScript) )
            // -- save the job id
            final jobId = (String)executor.parseJobId(result)
            updateStatus(jobId)
            log.debug "[${executor.name.toUpperCase()}] submitted process ${task.name} > jobId: $jobId; workDir: ${task.workDir}"

        }
        catch( Exception e ) {
            // update task exit status and message
            if( e instanceof ProcessNonZeroExitStatusException ) {
                task.exitStatus = e.getExitStatus()
                task.stdout = e.getReason()
                task.script = e.getCommand()
            }
            else {
                task.script = builder ? CmdLineHelper.toLine(builder.command()) : null
            }
            status = COMPLETED
            throw new ProcessFailedException("Error submitting process '${task.name}' for execution", e )
        }
    }

    private void updateStatus(String jobId) {
        if( task instanceof TaskArrayRun ) {
            for( int i=0; i<task.children.size(); i++ ) {
                final handler = task.children[i] as GridTaskHandler
                final arrayTaskId = ((TaskArrayExecutor)executor).getArrayTaskId(jobId, i)
                handler.updateStatus(arrayTaskId)
            }
        }
        else {
            this.jobId = jobId
            this.status = SUBMITTED
        }
    }

    private long startedMillis

    private long exitTimestampMillis0 = System.currentTimeMillis()

    private ExitStatusAwaiter exitAwaiter

    /**
     * When a process terminated save its exit status into the file defined by #exitFile
     *
     * @return The int value contained in the exit file or {@code null} if the file does not exist. When the
     * file contains an invalid number return {@code Integer#MAX_VALUE}
     */
    protected Integer readExitStatus() {

        /*
         * when the file does not exist return null, to force the monitor to continue to wait
         */
        BasicFileAttributes exitAttrs = null
        if( !exitFile || !(exitAttrs=FileHelper.readAttributes(exitFile)) || !exitAttrs.lastModifiedTime()?.toMillis() ) {
            if( log.isTraceEnabled() ) {
                if( !exitFile )
                    log.trace "JobId `$jobId` exit file is null"
                else
                    log.trace "JobId `$jobId` exit file: ${exitFile.toUriString()} - lastModified: ${exitAttrs?.lastModifiedTime()} - size: ${exitAttrs?.size()}"
            }

            // -- grace period: don't query the scheduler until enough time has elapsed since job submission
            def elapsed = System.currentTimeMillis() - startedMillis
            if( executor.queueInterval && elapsed < executor.queueInterval.toMillis() * 2.5 ) {
                exitAwaiter.reset()
                return null
            }

            // -- fetch the job status before returning a result
            // -- if the job is active, it is still running and the exit file absence is expected
            if( executor.checkActiveStatus(jobId, queue) ) {
                // make sure to reset the awaiter if the task is still active -- see #927
                exitAwaiter.reset()
                return null
            }

            // -- job is not active and exit file is missing: force NFS metadata refresh and delegate
            //    timeout tracking to the shared awaiter (mirrors ExitStatusAwaiter two-phase logic)
            String workDirList = null
            if( Global.session && FileHelper.workDirIsSharedFS ) {
                /*
                 * When the file is in a NFS folder, list the parent path to force refresh of NFS metadata
                 * before the awaiter re-reads the attributes, avoiding false-negative cache hits.
                 * http://stackoverflow.com/questions/3833127/alternative-to-file-exists-in-java
                 * http://superuser.com/questions/422061/how-to-determine-whether-a-directory-is-on-an-nfs-mounted-drive
                 */
                workDirList = FileHelper.listDirectory(task.workDir)
            }

            final result = exitAwaiter.read(exitFile)
            if( result == Integer.MAX_VALUE ) {
                def errMessage = []
                errMessage << "Failed to get exit status for process ${this} -- exitStatusReadTimeoutMillis: $exitStatusReadTimeoutMillis"
                errMessage << "Current queue status:"
                errMessage << executor.dumpQueueStatus()?.indent('> ')
                errMessage << "Content of workDir: ${task.workDir}"
                errMessage << workDirList?.indent('> ')
                log.debug errMessage.join('\n')
            }
            return result
        }

        /*
         * file is present: delegate content parsing and empty-file retry to the shared awaiter
         */
        return exitAwaiter.read(exitFile)
    }

    @Override
    boolean checkIfRunning() {

        if( isSubmitted() ) {
            if( isStarted() ) {
                status = RUNNING
                // use local timestamp because files are created on remote nodes which
                // may not have a synchronized clock
                startedMillis = System.currentTimeMillis()
                return true
            }
        }

        return false
    }

    private boolean isStarted() {

        BasicFileAttributes attr
        if( startFile && (attr=FileHelper.readAttributes(startFile)) && attr.lastModifiedTime()?.toMillis() > 0  )
            return true

        // check if the jobId is tracked in the queue status
        if( executor.checkStartedStatus(jobId, queue) )
            return true

        // to avoid unnecessary pressure on the file system check the existence of
        // the exit file on only on a time-periodic basis
        def now = System.currentTimeMillis()
        if( now - exitTimestampMillis0 > exitStatusReadTimeoutMillis ) {
            exitTimestampMillis0 = now
            // fix issue #268
            if( exitFile && (attr=FileHelper.readAttributes(exitFile)) && attr.lastModifiedTime()?.toMillis() > 0  )
                return true
        }

        return false
    }

    @Override
    boolean checkIfCompleted() {

        // verify the exit file exists
        Integer exit
        if( isRunning() && (exit = readExitStatus()) != null ) {
            // finalize the task
            task.exitStatus = exit
            task.stdout = outputFile
            task.stderr = errorFile
            status = COMPLETED
            return true
        }
        // sanity check
        else if( !passSanityCheck() ) {
            log.debug "Task sanity check failed > $task"
            task.stdout = outputFile
            task.stderr = errorFile
            status = COMPLETED
            return true
        }

        return false
    }

    protected boolean passSanityCheck() {
        Throttle.after(sanityCheckInterval, true) {
            if( isCompleted() ) {
                return true
            }
            if( task.workDir.exists() ) {
                return true
            }
            // if the task is not complete (ie submitted or running)
            // AND the work-dir does not exist ==> something is wrong
            task.error = new ProcessException("Task work directory is missing (!)")
            // sanity check does not pass
            return false
        }
    }

    @Override
    protected void killTask() {
        if( batch ) {
            batch.collect(executor, jobId)
        }
        else {
            executor.killTask(jobId)
        }
    }

    protected StringBuilder toStringBuilder( StringBuilder builder ) {
        builder << "jobId: $jobId; "

        super.toStringBuilder(builder)
        final exitAttrs = FileHelper.readAttributes(exitFile)

        builder << " started: " << (startedMillis ? startedMillis : '-') << ';'
        builder << " exited: " << (exitAttrs ? exitAttrs.lastModifiedTime() : '-') << '; '

        return builder
    }

    /**
     * @return An {@link nextflow.trace.TraceRecord} instance holding task runtime information
     */
    @Override
    TraceRecord getTraceRecord() {
        def trace = super.getTraceRecord()
        trace.put('native_id', jobId)
        return trace
    }

}
