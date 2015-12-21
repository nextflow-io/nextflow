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
import static nextflow.processor.TaskStatus.COMPLETED
import static nextflow.processor.TaskStatus.RUNNING
import static nextflow.processor.TaskStatus.SUBMITTED

import java.nio.file.Path

import groovy.util.logging.Slf4j
import nextflow.exception.ProcessFailedException
import nextflow.exception.ProcessSubmitException
import nextflow.file.FileHelper
import nextflow.processor.TaskHandler
import nextflow.processor.TaskRun
import nextflow.trace.TraceRecord
import nextflow.util.CmdLineHelper
import nextflow.util.Duration
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

    final static private READ_TIMEOUT = Duration.of('270sec') // 4.5 minutes

    BatchCleanup batch

    GridTaskHandler( TaskRun task, AbstractGridExecutor executor ) {
        super(task)

        this.executor = executor
        this.startFile = task.workDir.resolve(TaskRun.CMD_START)
        this.exitFile = task.workDir.resolve(TaskRun.CMD_EXIT)
        this.outputFile = task.workDir.resolve(TaskRun.CMD_OUTFILE)
        this.errorFile = task.workDir.resolve(TaskRun.CMD_ERRFILE)
        this.wrapperFile = task.workDir.resolve(TaskRun.CMD_RUN)
        final timeout = executor.session?.getExitReadTimeout(executor.name, READ_TIMEOUT) ?: READ_TIMEOUT
        this.exitStatusReadTimeoutMillis = timeout.toMillis()
        this.queue = task.config?.queue
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

    /*
     * {@inheritDocs}
     */
    @Override
    void submit() {
        log.debug "Launching process > ${task.name} -- work folder: ${task.workDir}"
        // -- create the wrapper script
        executor.createBashWrapperBuilder(task).build()

        // -- start the execution and notify the event to the monitor
        final builder = createProcessBuilder()

        def exitStatus = 0
        String result = null

        try {
            final process = builder.start()

            try {
                // -- forward the job launcher script to the command stdin if required
                if( executor.pipeLauncherScript() ) {
                    process.out << wrapperFile.text
                    process.out.close()
                }

                // -- wait the the process completes
                result = process.text
                exitStatus = process.waitFor()
                log.trace "submit ${task.name} > exit: $exitStatus\n$result\n"

                if( exitStatus ) {
                    throw new ProcessSubmitException("Failed to submit job to grid scheduler for execution")
                }

                // save the JobId in the
                this.jobId = executor.parseJobId(result)
                this.status = SUBMITTED
            }
            finally {
                // make sure to release all resources
                process.in.closeQuietly()
                process.out.closeQuietly()
                process.err.closeQuietly()
                process.destroy()
            }

        }
        catch( Exception e ) {
            task.exitStatus = exitStatus
            task.script = CmdLineHelper.toLine(builder.command())
            task.stdout = result
            status = COMPLETED
            throw new ProcessFailedException("Error submitting process '${task.name}' for execution", e )
        }

    }


    private long startedMillis

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
            workDirList = listDirectory(task.workDir)
        }

        /*
         * when the file does not exist return null, to force the monitor to continue to wait
         */
        def exitAttrs
        if( !exitFile || !(exitAttrs=exitFile.readAttributes()) || !exitAttrs.lastModifiedTime()?.toMillis() ) {
            // -- fetch the job status before return a result
            final active = executor.checkActiveStatus(jobId, queue)

            // --
            def elapsed = System.currentTimeMillis() - startedMillis
            if( elapsed < executor.queueInterval.toMillis() * 2.5 ) {
                return null
            }

            // -- if the job is active, this means that it is still running and thus the exit file cannot exist
            //    returns null to continue to wait
            if( active )
                return null

            // -- if the job is not active, something is going wrong
            //  * before returning an error code make (due to NFS latency) the file status could be in a incoherent state
            if( !exitTimestampMillis1 ) {
                log.debug "Exit file does not exist and the job is not running for task: $this -- Try to wait before kill it"
                exitTimestampMillis1 = System.currentTimeMillis()
            }

            def delta = System.currentTimeMillis() - exitTimestampMillis1
            if( delta < exitStatusReadTimeoutMillis ) {
                return null
            }

            def errMessage = []
            errMessage << "Failed to get exist status for process ${this} -- exitStatusReadTimeoutMillis: $exitStatusReadTimeoutMillis; delta: $delta"
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
                log.warn "Unable to parse process exit file: $exitFile -- bad value: '$status'"
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
                log.debug "File is returning an empty content $this -- Try to wait a while and .. pray."
                exitTimestampMillis2 = System.currentTimeMillis()
            }

            def delta = System.currentTimeMillis() - exitTimestampMillis2
            if( delta < exitStatusReadTimeoutMillis ) {
                return null
            }
            log.warn "Unable to read command status from: $exitFile after $delta ms"
        }

        return Integer.MAX_VALUE
    }

    @Override
    boolean checkIfRunning() {

        if( isSubmitted() ) {

            def startAttr
            if( startFile && (startAttr=startFile.readAttributes()) && startAttr.lastModifiedTime()?.toMillis() > 0) {
                status = RUNNING
                // use local timestamp because files are created on remote nodes which
                // may not have a synchronized clock
                startedMillis = System.currentTimeMillis()
                return true
            }

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

        return false
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
        final exitAttrs = exitFile.readAttributes()

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
        trace.native_id = jobId
        return trace
    }


    /*
     * When the file is in a NFS folder in order to avoid false negative
     * list the content of the parent path to force refresh of NFS metadata
     * http://stackoverflow.com/questions/3833127/alternative-to-file-exists-in-java
     * http://superuser.com/questions/422061/how-to-determine-whether-a-directory-is-on-an-nfs-mounted-drive
     */
    private String listDirectory(Path path) {

        String result = null
        Process process = null
        try {
            process = Runtime.runtime.exec("ls -la ${path}")
            def listStatus = process.waitFor()
            if( listStatus>0 ) {
                log.debug "Can't list folder: ${path} -- Exit status: $listStatus"
            }
            else {
                result = process.text
            }
        }
        catch( IOException e ) {
            log.debug "Can't list folder: $path -- Cause: ${e.message ?: e.toString()}"
        }
        finally {
            process?.destroy()
        }

        return result
    }
}
