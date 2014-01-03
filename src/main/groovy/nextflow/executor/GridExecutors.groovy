/*
 * Copyright (c) 2012, the authors.
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

import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.attribute.PosixFilePermissions

import groovy.transform.InheritConstructors
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.exception.InvalidExitException
import nextflow.processor.TaskConfig
import nextflow.processor.TaskHandler
import nextflow.processor.TaskMonitor
import nextflow.processor.TaskPollingMonitor
import nextflow.processor.TaskRun
import nextflow.util.CmdLineHelper
import org.apache.commons.io.IOUtils
/**
 * Generic task processor executing a task through a grid facility
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
abstract class AbstractGridExecutor extends AbstractExecutor {

    /**
     * Create a a queue holder for this executor
     * @return
     */
    def TaskMonitor createTaskMonitor() {
        final queueSize = session.getQueueSize(name, 50)
        final pollInterval = session.getPollIntervalMillis(name, 1_000)
        log.debug "Creating executor queue with size: $queueSize; poll-interval: $pollInterval"

        return new TaskPollingMonitor(session, queueSize, pollInterval)
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
        // set the input (when available)
        bash.input = task.stdin
        bash.scratch = taskConfig.scratch

        // set the environment
        bash.environment = task.processor.getProcessEnvironment()
        bash.environment.putAll( task.getInputEnvironment() )

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

    abstract List<String> killTaskCommand(def jobId);

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

    GridTaskHandler( TaskRun task, TaskConfig config, AbstractGridExecutor executor ) {
        super(task, config)
        this.executor = executor
        this.startFile = task.getCmdStartedFile()
        this.exitFile = task.getCmdExitFile()
        this.outputFile = task.getCmdOutputFile()
        this.wrapperFile = task.getCmdWrapperFile()
        this.exitStatusReadTimeoutMillis = executor.session?.getGridExitReadTimeoutMillis( executor.name ) ?: 90_000
    }

    /*
     * {@inheritDocs}
     */
    @Override
    void submit() {

        // -- log the qsub command
        def cli = executor.getSubmitCommandLine(task, wrapperFile)
        log.debug "sub command > '${cli}' -- process: ${task.name}"

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
                if( exitStatus ) {
                    new IllegalStateException("Grid submit command returned an error exit status: $exitStatus")
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
                throw new InvalidExitException("Error submitting task '${task.name}' for execution", e )
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


    private long exitTimestampMillis

    private int exitEmptyCount

    /**
     * When a process terminated save its exit status into the file defined by #exitFile
     *
     * @return The int value contained in the exit file or {@code null} if the file does not exist. When the
     * file contains an invalid number return {@code Integer#MAX_VALUE}
     */
    private Integer readExitStatus() {

        /*
         * when the file does not exist return null, to force the monitor to continue to wait
         */
        if( !exitFile || !exitFile.exists() ) {
            return null
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
            if( !exitTimestampMillis ) {
                exitTimestampMillis = System.currentTimeMillis()
            }

            def delta = System.currentTimeMillis() - exitTimestampMillis
            if( delta < exitStatusReadTimeoutMillis ) {
                log.debug "File is returning an empty content ($exitEmptyCount): $exitFile -- Try to wait a while and .. pray."
                return null
            }
            log.warn "Unable to read command status from: $exitFile after $delta ms"
        }

        return Integer.MAX_VALUE
    }

    @Override
    boolean checkIfRunning() {

        if( isSubmitted() && startFile && startFile.exists() ) {
            status = RUNNING
            return true
        }

        return false
    }

    @Override
    boolean checkIfCompleted() {

        def _exit
        if( isRunning() && (_exit = readExitStatus()) != null ) {
            // finalize the task
            task.exitStatus = _exit
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

        // the requested queue name
        if( taskConfig.queue ) {
            result << '-q'  << taskConfig.queue
        }

        // max task duration
        if( taskConfig.maxDuration ) {
            result << '-l' << "h_rt=${taskConfig.maxDuration.format('HH:mm:ss')}"
        }

        // task max memory
        if( taskConfig.maxMemory ) {
            result << '-l' << "virtual_free=${taskConfig.maxMemory.toString().replaceAll(/[\sB]/,'')}"
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
        return lines ? lines[-1].trim() : null
    }


    @PackageScope
    List<String> killTaskCommand(jobId) {
        ['qdel', '-j', jobId?.toString()]
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
        Files.setPosixFilePermissions( scriptFile, PosixFilePermissions.fromString('rwx------'))

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

        new IllegalStateException("Not a valid 'bsub' output:\n$text");
    }

    @Override
    List<String> killTaskCommand( def jobId ) {
        ['bkill', jobId?.toString() ]
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

        throw new IllegalStateException()
    }


    List<String> killTaskCommand(def jobId) {
        ['scancel', jobId?.toString() ]
    }
}