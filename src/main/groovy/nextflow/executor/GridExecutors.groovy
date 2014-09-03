/*
 * Copyright (c) 2013-2014, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2014, Paolo Di Tommaso and the respective authors.
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
import static nextflow.processor.TaskHandler.Status.COMPLETED
import static nextflow.processor.TaskHandler.Status.RUNNING
import static nextflow.processor.TaskHandler.Status.SUBMITTED

import java.nio.file.Path

import groovy.transform.InheritConstructors
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.exception.ProcessFailedException
import nextflow.exception.ProcessSubmitException
import nextflow.processor.TaskConfig
import nextflow.processor.TaskHandler
import nextflow.processor.TaskMonitor
import nextflow.processor.TaskPollingMonitor
import nextflow.processor.TaskRun
import nextflow.util.CmdLineHelper
import nextflow.util.Duration
import org.apache.commons.io.IOUtils
import org.apache.commons.lang.StringUtils
/**
 * Generic task processor executing a task through a grid facility
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
abstract class AbstractGridExecutor extends Executor {


    protected Duration queueInterval

    /**
     * Initialize the executor class
     */
    void init() {
        super.init()
        queueInterval = session.getQueueStatInterval(name)
        log.debug "Creating executor '$name' > queue-stat-interval: ${queueInterval}"
    }

    /**
     * Create a a queue holder for this executor
     * @return
     */
    def TaskMonitor createTaskMonitor() {
        return TaskPollingMonitor.create(session, name, 50, Duration.of('1 sec'))
    }

    /*
     * Prepare and launch the task in the underlying execution platform
     */
    def GridTaskHandler createTaskHandler(TaskRun task) {
        assert task
        assert task.workDirectory

        final folder = task.workDirectory
        log.debug "Launching process > ${task.name} -- work folder: $folder"

        final bash = new BashWrapperBuilder(task)

        // staging input files
        bash.stagingScript = {
            final files = task.getInputFiles()
            final staging = stagingFilesScript(files)
            return staging
        }

        // unstage script
        bash.unstagingScript = {
            return unstageOutputFilesScript(task)
        }

        // create the wrapper script
        bash.build()

        return new GridTaskHandler(task, taskConfig, this)
    }


    /**
     * Build up the platform native command line used to submit the job wrapper
     * execution request to the underlying grid, e.g. {@code qsub -q something script.job}
     *
     * @param task The task instance descriptor
     * @return A list holding the command line
     */
    abstract List<String> getSubmitCommandLine(TaskRun task, Path scriptFile )

    /**
     * Given the string returned the by grid submit command, extract the process handle i.e. the grid jobId
     */
    abstract parseJobId( String text );

    /**
     * Kill a grid job
     *
     * @param jobId The ID of the job to kill
     */
    final void killTask( def jobId )  {
        new ProcessBuilder(killTaskCommand(jobId)).start()
    }

    /**
     * The command to be used to kill a grid job
     * @param jobId The job ID to be kill
     * @return The command line to be used to kill the specified job
     */
    protected abstract List<String> killTaskCommand(def jobId);

    /**
     * @return Parse the {@code clusterOptions} configuration option and return the entries as a list of values
     */
    final protected List<String> getClusterOptionsAsList() {

        if ( !taskConfig.clusterOptions ) {
            return null
        }

        if( taskConfig.clusterOptions instanceof Collection ) {
            return new ArrayList<String>(taskConfig.clusterOptions as Collection)
        }
        else {
            return CmdLineHelper.splitter( taskConfig.clusterOptions.toString() )
        }
    }

    /**
     * @return Parse the {@code clusterOptions} configuration option and return the entries as a string
     */
    final protected String getClusterOptionsAsString() {

        if( !taskConfig.clusterOptions ) {
            return null
        }

        def value = taskConfig.clusterOptions
        value instanceof Collection ? value.join(' ') : value.toString()
    }

    /**
     * Status as returned by the grid engine
     */
    static protected enum QueueStatus { PENDING, RUNNING, HOLD, ERROR, DONE, UNKNWON }

    /**
     * @return The status for all the scheduled and running jobs
     */
    Map<?,QueueStatus> getQueueStatus() {

        List cmd = queueStatusCommand( taskConfig.queue )
        if( !cmd ) return null

        try {
            log.trace "Getting grid queue status: ${cmd.join(' ')}"

            def process = new ProcessBuilder(cmd).start()
            def result = process.text
            process.waitForOrKill( 10 * 1000 )
            def exit = process.exitValue()

            log.trace "status result > exit: $exit\n$result\n"

            return ( exit == 0 ) ? parseQueueStatus( result ) : null

        }
        catch( Exception e ) {
            log.warn "Unable to fetch queue status -- See the log file for details", e
            return null
        }

    }

    @PackageScope
    String dumpQueueStatus() {
        def result = new StringBuilder()
        fQueueStatus?.each { k, v ->
            result << '  job: ' << StringUtils.leftPad(k?.toString(),6) << ': ' << v?.toString() << '\n'
        }
        return result.toString()
    }

    /**
     * @param queue The command for which the status of jobs has be to read
     * @return The command line to be used to retried the job statuses
     */
    protected abstract List<String> queueStatusCommand(queue)

    /**
     * Parse the command stdout produced by the command line {@code #queueStatusCommand}
     * @param text
     * @return
     */
    protected abstract Map<?,QueueStatus> parseQueueStatus( String text )

    /**
     * Store jobs status
     */
    protected Map<Object,QueueStatus> fQueueStatus = null

    /**
     * Verify that a job in a 'active' state i.e. RUNNING or HOLD
     *
     * @param jobId The job for which verify the status
     * @return {@code true} if the job is in RUNNING or HOLD status, or even if it is temporarily unable
     *  to retrieve the job status for some
     */
    public boolean checkActiveStatus( jobId ) {

        // -- fetch the queue status
        fQueueStatus = queueInterval.throttle(null) { getQueueStatus() }
        if( fQueueStatus == null ) // no data is returned, so return true
            return true

        if( !fQueueStatus.containsKey(jobId) )
            return false

        return fQueueStatus[jobId] == QueueStatus.RUNNING || fQueueStatus[jobId] == QueueStatus.HOLD

    }

}


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

    /** The wrapper file used to execute the user script */
    final Path wrapperFile

    /** The unique job ID as provided by the underlying grid platform */
    private jobId

    private long exitStatusReadTimeoutMillis

    final static private READ_TIMEOUT = Duration.of('270sec') // 4.5 minutes

    GridTaskHandler( TaskRun task, TaskConfig config, AbstractGridExecutor executor ) {
        super(task, config)

        this.executor = executor
        this.startFile = task.getCmdStartedFile()
        this.exitFile = task.getCmdExitFile()
        this.outputFile = task.getCmdOutputFile()
        this.wrapperFile = task.getCmdWrapperFile()
        final timeout = executor.session?.getExitReadTimeout(executor.name, READ_TIMEOUT) ?: READ_TIMEOUT
        this.exitStatusReadTimeoutMillis = timeout.toMillis()
    }

    /*
     * {@inheritDocs}
     */
    @Override
    void submit() {

        // -- log the qsub command
        def cli = executor.getSubmitCommandLine(task, wrapperFile)
        log.trace "submit ${task.name} > cli: ${cli}"

        /*
         * launch 'sub' script wrapper
         */
        ProcessBuilder builder = new ProcessBuilder()
                .directory(task.workDirectory.toFile())
                .command( cli as String[] )
                .redirectErrorStream(true)

        // -- start the execution and notify the event to the monitor
        Process process = builder.start()

        try {
            def exitStatus = 0
            String result = null
            try {
                // -- wait the the process completes
                result = process.text
                exitStatus = process.waitFor()
                log.trace "submit ${task.name} > exit: $exitStatus\n$result\n"

                if( exitStatus ) {
                    throw new ProcessSubmitException("Failed to submit job to the grid engine job for execution")
                }

                // save the JobId in the
                this.jobId = executor.parseJobId(result)
                this.status = SUBMITTED
            }
            catch( Exception e ) {
                task.exitStatus = exitStatus
                task.script = CmdLineHelper.toLine(cli)
                task.stdout = result
                status = COMPLETED
                throw new ProcessFailedException("Error submitting process '${task.name}' for execution", e )
            }

        }
        finally {
            // make sure to release all resources
            IOUtils.closeQuietly(process.in)
            IOUtils.closeQuietly(process.out)
            IOUtils.closeQuietly(process.err)
            process.destroy()
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

        /*
         * when the file does not exist return null, to force the monitor to continue to wait
         */
        if( !exitFile || !exitFile.exists() || !exitFile.lastModified() ) {
            // -- fetch the job status before return a result
            final active = executor.checkActiveStatus(jobId)

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

            log.debug "Failed to get exist status for process ${this} -- exitStatusReadTimeoutMillis: $exitStatusReadTimeoutMillis; delta: $delta"

            // -- dump current queue stats
            log.debug "Current queue status:\n" + executor.dumpQueueStatus()

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

        if( isSubmitted()  ) {

            if( startFile && startFile.exists() && startFile.lastModified() > 0) {
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
            status = COMPLETED
            return true

        }

        return false

    }

    @Override
    void kill() {
        executor.killTask(jobId)
    }

    protected StringBuilder toStringBuilder( StringBuilder builder ) {
        builder << "jobId: $jobId; "

        super.toStringBuilder(builder)

        builder << " started: " << (startedMillis ? startedMillis : '-') << ';'
        builder << " exited: " << (exitFile.exists() ? exitFile.lastModified() : '-') << '; '

        return builder
    }

}

/**
 * Execute a task script by running it on the SGE/OGE cluster
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class SgeExecutor extends AbstractGridExecutor {


    /*
     * Prepare the 'qsub' cmdline. The following options are used
     * - wd: define the job working directory
     * - terse: output just the job id on the output stream
     * - j: redirect qsub stdout and stderr to the same file (join)
     * - sync: wait for the job completion
     * -V: export the current environment
     */
    List<String> getSubmitCommandLine(TaskRun task, Path scriptFile ) {

        final result = new ArrayList<String>()

        result << 'qsub'
        result << '-wd' << task.workDirectory?.toString()
        result << '-N' << "nf-${task.name.replace(' ','_')}"
        result << '-o' << '/dev/null'
        result << '-j' << 'y'
        result << '-terse'
        result << '-V'
        /*
         * By using command line option -notify SIGUSR1 will be sent to your script prior to SIGSTOP
         * and SIGUSR2 will be sent to your script prior to SIGKILL
         */
        result << '-notify'

        // the requested queue name
        if( taskConfig.queue ) {
            result << '-q'  << taskConfig.queue
        }

        // max task duration
        if( taskConfig.maxDuration ) {
            final duration = taskConfig.maxDuration as Duration
            result << "-l" << "h_rt=${duration.format('HH:mm:ss')}"
        }

        // task max memory
        if( taskConfig.maxMemory ) {
            result << "-l" << "virtual_free=${taskConfig.maxMemory.toString().replaceAll(/[\sB]/,'')}"
        }

        // -- at the end append the command script wrapped file name
        if( taskConfig.clusterOptions ) {
            result.addAll( getClusterOptionsAsList() )
        }

        // -- last entry to 'script' file name
        result << scriptFile.getName()

        return result
    }

    @Override
    def parseJobId( String text ) {
        // return always the last line
        def lines = text.trim().readLines()
        def id = lines[-1].trim()
        if( id?.toString()?.isInteger() )
            return id

        throw new IllegalStateException("Invalid SGE submit response:\n$text\n\n")
    }


    @PackageScope
    List<String> killTaskCommand(jobId) {
        ['qdel', '-j', jobId?.toString()]
    }

    @Override
    protected List<String> queueStatusCommand(Object queue) {
        def result = ['qstat']
        if( queue )
            result << '-q' << queue.toString()

        return result
    }

    static private Map DECODE_STATUS = [
            'r': QueueStatus.RUNNING,
            'qw': QueueStatus.PENDING,
            'hqw': QueueStatus.HOLD,
            'Eqw': QueueStatus.ERROR
    ]

    @Override
    protected Map<?, QueueStatus> parseQueueStatus(String text) {

        def result = [:]
        text?.eachLine{ String row, int index ->
            if( index< 2 ) return
            def cols = row.split(/\s+/)
            if( cols.size()>5 ) {
                result.put( cols[0], DECODE_STATUS[cols[4]] )
            }
        }

        return result
    }

}




/**
 * Processor for LSF resource manager (DRAFT)
 *
 * See http://en.wikipedia.org/wiki/Platform_LSF
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class LsfExecutor extends AbstractGridExecutor {

    @Override
    List<String> getSubmitCommandLine(TaskRun task, Path scriptFile ) {

        // note: LSF requires the job script file to be executable
        scriptFile.setPermissions(7,0,0)

        final result = new ArrayList<String>()
        result << 'bsub'
        result << '-cwd' << task.workDirectory?.toString()
        result << '-o' << '/dev/null'

        // add other parameters (if any)
        if( taskConfig.queue ) {
            result << '-q'  << taskConfig.queue
        }

        // -- the job name
        result << '-J' << "nf-${task.name.replace(' ','_')}"

        // -- at the end append the command script wrapped file name
        if( taskConfig.clusterOptions ) {
            result.addAll( getClusterOptionsAsList() )
        }

        // -- last entry to 'script' file name
        result << "./${scriptFile.getName()}"

        return result

    }

    @Override
    def parseJobId(String text) {

        def pattern = ~/Job <(\d+)> is submitted/
        for( String line : text.readLines() ) {
            def m = pattern.matcher(line)
            if( m.find() ) {
                return m[0][1].toString()
            }
        }

        new IllegalStateException("Invalid LSF submit response:\n$text\n\n");
    }

    @Override
    List<String> killTaskCommand( def jobId ) {
        ['bkill', jobId?.toString() ]
    }

    @Override
    protected List<String> queueStatusCommand(Object queue) {

        def result = ['bjobs', '-o',  'JOBID STAT SUBMIT_TIME delimiter=\',\'', '-noheader']

        if( queue )
            result << '-q' << queue

        return result

    }

    private static Map DECODE_STATUS = [
        'PEND': QueueStatus.PENDING,
        'RUN': QueueStatus.RUNNING,
        'PSUSP': QueueStatus.HOLD,
        'USUSP': QueueStatus.HOLD,
        'SSUSP': QueueStatus.HOLD,
        'DONE': QueueStatus.DONE,
        'EXIT': QueueStatus.ERROR,
        'UNKWN': QueueStatus.ERROR,
        'ZOMBI': QueueStatus.ERROR,
    ]

    @Override
    protected Map<?, QueueStatus> parseQueueStatus(String text) {

        def result = [:]

        text.eachLine { String line ->
            def cols = line.split(',')
            if( cols.size() == 3 ) {
                result.put( cols[0], DECODE_STATUS.get(cols[1]) )
            }
        }

        return result
    }
}

/**
 * Processor for SLURM resource manager (DRAFT)
 *
 * See http://computing.llnl.gov/linux/slurm/
 *
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@InheritConstructors
class SlurmExecutor extends AbstractGridExecutor {


    /**
     * -c, --cpus-per-task=ncpus   number of cpus required per task
     * -D, --chdir=path            change remote current working directory
     * -e, --error=err             location of stderr redirection
     * -E, --preserve-env          env vars for node and task counts override
     * -i, --input=in              location of stdin redirection
     * -J, --job-name=jobname      name of job
     * -o, --output=out            location of stdout redirection
     * -Q, --quiet                 quiet mode (suppress informational messages)
     * -t, --time=minutes          time limit
     * -u, --unbuffered            do not line-buffer stdout/err
     * --mem=MB                minimum amount of real memory
     * --mincpus=n             minimum number of logical processors (threads) per node
     * --tmp=MB                minimum amount of temporary disk
     * --mem-per-cpu=MB        maximum amount of real memory per allocated cpu required by the job. --mem >= --mem-per-cpu if --mem is specified.
     *
     * @param task
     * @return
     */
    @Override
    List<String> getSubmitCommandLine(TaskRun task, Path scriptFile ) {

        final result = new ArrayList<String>()

        result << 'sbatch'
        result << '-D' << task.workDirectory.toString()
        result << '-J' << "nf-${task.name.replace(' ','_')}"
        result << '-o' << '/dev/null'

        if( taskConfig.maxDuration ) {
            result << '-t' << taskConfig.maxDuration.format('HH:mm:ss')
        }

        // -- at the end append the command script wrapped file name
        if( taskConfig.clusterOptions ) {
            result.addAll( getClusterOptionsAsList() )
        }

        // -- last entry to 'script' file name
        // replace with the 'shell' attribute
        result << scriptFile.getName()

    }

    @Override
    def parseJobId(String text) {
        def pattern = ~ /Submitted batch job (\d+)/
        for( String line : text.readLines() ) {
            def m = pattern.matcher(line)
            if( m.matches() ) {
                return m[0][1].toString()
            }
        }

        throw new IllegalStateException("Invalid SLURM submit response:\n$text\n\n")
    }


    protected List<String> killTaskCommand(def jobId) {
        ['scancel', jobId?.toString() ]
    }

    @Override
    protected List<String> queueStatusCommand(Object queue) {
        if( queue )
            log.debug "SLURM executor does not support queue parameter on queue status"

        return ['squeue','-h','-o \'%i %t\'']
    }

    static private Map STATUS_MAP = [
            'PD': QueueStatus.PENDING,  // (pending)
            'R': QueueStatus.RUNNING,   // (running)
            'CA': QueueStatus.ERROR,    // (cancelled)
            'CF': QueueStatus.PENDING,  // (configuring)
            'CG': QueueStatus.RUNNING,  // (completing)
            'CD': QueueStatus.DONE,     // (completed)
            'F': QueueStatus.ERROR,     // (failed),
            'TO': QueueStatus.ERROR,    // (timeout),
            'NF': QueueStatus.ERROR     // (node failure)
    ]

    @Override
    protected Map<?, QueueStatus> parseQueueStatus(String text) {

        def result = [:]

        text.eachLine { String line ->
            def cols = line.split(/\s+/)
            if( cols.size() == 2 ) {
                result.put( cols[0], STATUS_MAP.get(cols[1]) )
            }
        }

        return result
    }
}


@Slf4j
@InheritConstructors
class PbsExecutor extends AbstractGridExecutor {

    List<String> getSubmitCommandLine(TaskRun task, Path scriptFile ) {

        final result = new ArrayList<String>()

        result << 'qsub'
        result << '-d' << task.workDirectory?.toString()
        result << '-N' << "nf-${task.name.replace(' ','_')}"
        result << '-o' << '/dev/null'
        result << '-e' << '/dev/null'
        result << '-V'

        // the requested queue name
        if( taskConfig.queue ) {
            result << '-q'  << (String)taskConfig.queue
        }

//        // max task duration
//        if( taskConfig.maxDuration ) {
//            final duration = taskConfig.maxDuration as Duration
//            result << "-l" << "h_rt=${duration.format('HH:mm:ss')}"
//        }
//
        // task max memory
        if( taskConfig.maxMemory ) {
            result << "-l" << "mem=${taskConfig.maxMemory.toString().replaceAll(/[\s]/,'')}"
        }

        // -- at the end append the command script wrapped file name
        if( taskConfig.clusterOptions ) {
            result.addAll( getClusterOptionsAsList() )
        }

        // -- last entry to 'script' file name
        result << scriptFile.getName()

        return result
    }

    @Override
    def parseJobId( String text ) {
        // return always the last line
        def result = text?.trim()
        if( result )
            return result

        throw new IllegalStateException("Invalid PBS/Torque submit response:\n$text\n\n")
    }


    @PackageScope
    List<String> killTaskCommand(jobId) {
        ['qdel', jobId?.toString()]
    }

    @Override
    protected List<String> queueStatusCommand(Object queue) {
        def result = ['qstat']
        if( queue )
            result << queue.toString()

        return result
    }

    static private Map DECODE_STATUS = [
            'C': QueueStatus.DONE,
            'R': QueueStatus.RUNNING,
            'Q': QueueStatus.PENDING,
            'H': QueueStatus.HOLD,
            'S': QueueStatus.HOLD
    ]

    @Override
    protected Map<?, QueueStatus> parseQueueStatus(String text) {

        def result = [:]
        text?.eachLine{ String row, int index ->
            if( index< 2 ) return
            def cols = row.split(/\s+/)
            if( cols.size()>5 ) {
                result.put( cols[0], DECODE_STATUS[cols[4]] ?: AbstractGridExecutor.QueueStatus.UNKNWON )
            }
        }

        return result
    }


}
