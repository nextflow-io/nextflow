/*
 * Copyright 2020-2022, Seqera Labs
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

import static nextflow.processor.TaskStatus.*

import java.nio.file.Path
import java.nio.file.attribute.BasicFileAttributes
import java.time.temporal.ChronoUnit
import java.util.function.Predicate
import java.util.regex.Pattern

import dev.failsafe.Failsafe
import dev.failsafe.RetryPolicy
import dev.failsafe.event.EventListener
import dev.failsafe.event.ExecutionAttemptedEvent
import dev.failsafe.function.CheckedSupplier
import groovy.transform.Memoized
import groovy.util.logging.Slf4j
import nextflow.exception.ProcessNonZeroExitStatusException
import nextflow.exception.ProcessException
import nextflow.exception.ProcessFailedException
import nextflow.file.FileHelper
import nextflow.processor.TaskHandler
import nextflow.processor.TaskRun
import nextflow.trace.TraceRecord
import nextflow.util.CmdLineHelper
import nextflow.util.Duration
import nextflow.util.Throttle
/**
 * Handles a job execution in the underlying grid platform
 */
@Slf4j
class GridTaskHandler extends TaskHandler {

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

    final static private READ_TIMEOUT = Duration.of('270sec') // 4.5 minutes

    BatchCleanup batch

    /** only for testing purpose */
    protected GridTaskHandler() {}

    GridTaskHandler( TaskRun task, AbstractGridExecutor executor ) {
        super(task)

        this.executor = executor
        this.startFile = task.workDir.resolve(TaskRun.CMD_START)
        this.exitFile = task.workDir.resolve(TaskRun.CMD_EXIT)
        this.outputFile = task.workDir.resolve(TaskRun.CMD_OUTFILE)
        this.errorFile = task.workDir.resolve(TaskRun.CMD_ERRFILE)
        this.wrapperFile = task.workDir.resolve(TaskRun.CMD_RUN)
        final duration = executor.session?.getExitReadTimeout(executor.name, READ_TIMEOUT) ?: READ_TIMEOUT
        this.exitStatusReadTimeoutMillis = duration.toMillis()
        this.queue = task.config?.queue
        this.sanityCheckInterval = duration
    }

    protected ProcessBuilder createProcessBuilder() {

        // -- log the qsub command
        def cli = executor.getSubmitCommandLine(task, wrapperFile)
        log.trace "start process ${task.name} > cli: ${cli}"

        /*
         * launch 'sub' script wrapper
         */
        ProcessBuilder builder = new ProcessBuilder()
                .command( cli as String[] )
                .redirectErrorStream(true)
                .directory(task.workDir.toFile())

        return builder
    }

    @Memoized
    protected Predicate<? extends Throwable> retryCondition(String reasonPattern) {
        final pattern = Pattern.compile(reasonPattern)
        return new Predicate<Throwable>() {
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

        final delay = executor.session.getConfigAttribute("executor.retry.delay", '500ms') as Duration
        final maxDelay = executor.session.getConfigAttribute("executor.retry.maxDelay", '30s') as Duration
        final jitter = executor.session.getConfigAttribute("executor.retry.jitter", '0.25') as double
        final maxAttempts = executor.session.getConfigAttribute("executor.retry.maxAttempts", '3') as int
        final reason = executor.session.getConfigAttribute("executor.submit.retry.reason", 'Socket timed out') as String

        final listener = new EventListener<ExecutionAttemptedEvent>() {
            @Override
            void accept(ExecutionAttemptedEvent event) throws Throwable {
                final failure = event.getLastFailure()
                if( failure instanceof ProcessNonZeroExitStatusException ) {
                    final msg = """\
                        Failed to submit process '${task.name}'
                         - attempt : ${event.attemptCount}
                         - command : ${CmdLineHelper.toLine(failure.command)}
                         - reason  : $failure.reason
                        """.stripIndent()
                    log.warn msg

                } else {
                    log.debug("Unexpected retry failure: ${failure?.message}", failure)
                }
            }
        }

        return RetryPolicy.<T>builder()
                .handleIf(retryCondition(reason))
                .withBackoff(delay.toMillis(), maxDelay.toMillis(), ChronoUnit.MILLIS)
                .withMaxAttempts(maxAttempts)
                .withJitter(jitter)
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
                process.out << pipeScript
                process.out.close()
            }

            // -- wait the the process completes
            final result = process.text
            final exitStatus = process.waitFor()

            if( exitStatus ) {
                throw new ProcessNonZeroExitStatusException("Failed to submit process to grid scheduler for execution", result, exitStatus, builder.command())
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
    
    /*
     * {@inheritDocs}
     */
    @Override
    void submit() {
        ProcessBuilder builder = null
        try {
            // -- create the wrapper script
            executor.createBashWrapperBuilder(task).build()
            // -- start the execution and notify the event to the monitor
            builder = createProcessBuilder()
            // -- forward the job launcher script to the command stdin if required
            final stdinScript = executor.pipeLauncherScript() ? wrapperFile.text : null
            // -- execute with a re-triable strategy
            final result = safeExecute( () -> processStart(builder, stdinScript) )
            // -- save the JobId in the
            this.jobId = executor.parseJobId(result)
            this.status = SUBMITTED
            log.debug "[${executor.name.toUpperCase()}] submitted process ${task.name} > jobId: $jobId; workDir: ${task.workDir}"

        }
        catch( Exception e ) {
            // update task exit status and message
            if( e instanceof ProcessNonZeroExitStatusException ) {
                task.exitStatus = e.getExitStatus()
                task.stdout = e.getReason()
            }
            task.script = builder ? CmdLineHelper.toLine(builder.command()) : null
            status = COMPLETED
            throw new ProcessFailedException("Error submitting process '${task.name}' for execution", e )
        }

    }


    private long startedMillis

    private long exitTimestampMillis0 = System.currentTimeMillis()

    private long exitTimestampMillis1

    private long exitTimestampMillis2

    /**
     * When a process terminated save its exit status into the file defined by #exitFile
     *
     * @return The int value contained in the exit file or {@code null} if the file does not exist. When the
     * file contains an invalid number return {@code Integer#MAX_VALUE}
     */
    protected Integer readExitStatus() {

        String workDirList = null
        if( exitTimestampMillis1 && FileHelper.workDirIsNFS ) {
            /*
             * When the file is in a NFS folder in order to avoid false negative
             * list the content of the parent path to force refresh of NFS metadata
             * http://stackoverflow.com/questions/3833127/alternative-to-file-exists-in-java
             * http://superuser.com/questions/422061/how-to-determine-whether-a-directory-is-on-an-nfs-mounted-drive
             */
            workDirList = FileHelper.listDirectory(task.workDir)
        }

        /*
         * when the file does not exist return null, to force the monitor to continue to wait
         */
        def exitAttrs = null
        if( !exitFile || !(exitAttrs=FileHelper.readAttributes(exitFile)) || !exitAttrs.lastModifiedTime()?.toMillis() ) {
            if( log.isTraceEnabled() ) {
                if( !exitFile )
                    log.trace "JobId `$jobId` exit file is null"
                else
                    log.trace "JobId `$jobId` exit file: ${exitFile.toUriString()} - lastModified: ${exitAttrs?.lastModifiedTime()} - size: ${exitAttrs?.size()}"
            }
            // -- fetch the job status before return a result
            final active = executor.checkActiveStatus(jobId, queue)

            // --
            def elapsed = System.currentTimeMillis() - startedMillis
            if( elapsed < executor.queueInterval.toMillis() * 2.5 ) {
                return null
            }

            // -- if the job is active, this means that it is still running and thus the exit file cannot exist
            //    returns null to continue to wait
            if( active ) {
                // make sure to reset exit time if the task is active -- see #927
                exitTimestampMillis1 = 0
                return null
            }

            // -- if the job is not active, something is going wrong
            //  * before returning an error code make (due to NFS latency) the file status could be in a incoherent state
            if( !exitTimestampMillis1 ) {
                log.trace "Exit file does not exist for and the job is not running for task: $this -- Try to wait before kill it"
                exitTimestampMillis1 = System.currentTimeMillis()
            }

            def delta = System.currentTimeMillis() - exitTimestampMillis1
            if( delta < exitStatusReadTimeoutMillis ) {
                return null
            }

            def errMessage = []
            errMessage << "Failed to get exit status for process ${this} -- exitStatusReadTimeoutMillis: $exitStatusReadTimeoutMillis; delta: $delta"
            // -- dump current queue stats
            errMessage << "Current queue status:"
            errMessage << executor.dumpQueueStatus()?.indent('> ')
            // -- dump directory listing
            errMessage << "Content of workDir: ${task.workDir}"
            errMessage << workDirList?.indent('> ')
            log.debug errMessage.join('\n')

            return Integer.MAX_VALUE
        }

        /*
         * read the exit file, it should contain the executed process exit status
         */
        def status = exitFile.text?.trim()
        if( status ) {
            try {
                return status.toInteger()
            }
            catch( Exception e ) {
                log.warn "Unable to parse process exit file: ${exitFile.toUriString()} -- bad value: '$status'"
            }
        }

        else {
            /*
             * Since working with NFS it may happen that the file exists BUT it is empty due to network latencies,
             * before retuning an invalid exit code, wait some seconds.
             *
             * More in detail:
             * 1) the very first time that arrive here initialize the 'exitTimestampMillis' to the current timestamp
             * 2) when the file is empty but less than 5 seconds are spent from the first check, return null
             *    this will force the monitor to continue to wait for job termination
             * 3) if more than 5 seconds are spent, and the file is empty return MAX_INT as an invalid exit status
             *
             */
            if( !exitTimestampMillis2 ) {
                log.debug "File is returning empty content: $this -- Try to wait a while... and pray."
                exitTimestampMillis2 = System.currentTimeMillis()
            }

            def delta = System.currentTimeMillis() - exitTimestampMillis2
            if( delta < exitStatusReadTimeoutMillis ) {
                return null
            }
            log.warn "Unable to read command status from: ${exitFile.toUriString()} after $delta ms"
        }

        return Integer.MAX_VALUE
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
        def exit
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
            // AND the work-dir does not exists ==> something is wrong
            task.error = new ProcessException("Task work directory is missing (!)")
            // sanity check does not pass
            return false
        }
    }

    @Override
    void kill() {
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
    public TraceRecord getTraceRecord() {
        def trace = super.getTraceRecord()
        trace.put('native_id', jobId)
        return trace
    }

}
