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
import groovy.util.logging.Slf4j
import nextflow.processor.TaskRun
import nextflow.util.ByteDumper
import nextflow.util.CmdLineHelper
import nextflow.util.Duration
import org.apache.commons.io.IOUtils

/**
 * Generic task processor executing a task through a grid facility
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
abstract class AbstractGridExecutor extends AbstractExecutor {

    protected static final COMMAND_SCRIPT_FILENAME = '.command.sh'

    protected static final COMMAND_OUTPUT_FILENAME = '.command.out'

    protected static final COMMAND_INPUT_FILE = '.command.in'

    protected static final COMMAND_ENV_FILENAME = '.command.env'

    protected static final JOB_OUT_FILENAME = '.job.out'

    protected static final JOB_SCRIPT_FILENAME = '.job.run'


    /*
     * Prepare and launch the task in the underlying execution platform
     */
    @Override
    void launchTask( TaskRun task ) {
        assert task
        assert task.workDirectory

        final folder = task.workDirectory
        log.debug "Lauching task > ${task.name} -- work folder: $folder"

        /*
         * save the environment to a file
         */
        def envFile = new File(folder, COMMAND_ENV_FILENAME)
        createEnvironmentFile(task, envFile)

        /*
         * save the command input (if any)
         * the file content will be piped to the executed user script
         */
        File cmdInputFile = createCommandInputFile(task)

        /*
         * save the 'user' script to be executed
         */
        def scriptFile = createCommandScriptFile(task)


        /*
         * create the job wrapper script file
         */
        def cmdOutFile = new File(folder, COMMAND_OUTPUT_FILENAME)
        def runnerFile = createJobWrapperFile(task, scriptFile, envFile, cmdInputFile, cmdOutFile)


        /*
         * Finally submit the job script for execution
         */
        submitJob(task, runnerFile, cmdOutFile)

    }

    /**
     * Save the main task script to a file having executable permission
     */
    protected File createCommandScriptFile( TaskRun task ) {
        assert task

        def scriptFile = new File(task.workDirectory, COMMAND_SCRIPT_FILENAME)
        scriptFile.text = task.processor.normalizeScript(task.script.toString())
        scriptFile.setExecutable(true)
        task.script = scriptFile

        return scriptFile
    }

    /**
     * Save the task {@code stdin} to a file to be piped when it is executed
     *
     * @param task
     * @return
     */
    protected File createCommandInputFile( TaskRun task ) {

        if( task.input == null ) {
            return null
        }

        def result = new File( task.workDirectory, COMMAND_INPUT_FILE )
        result.text = task.input
        return result
    }

    /**
     * Create a script which wraps the target script to be able to manage custom environment variable, scratch directory, etc
     *
     * @param task The task instance descriptor
     * @param scriptFile The task main script file
     * @param envFile The file containing the task environment
     * @param cmdInputFile The 'stdin' file (if any)
     * @param cmdOutFile The file that where save the task output
     * @return The script wrapper file
     */
    protected File createJobWrapperFile( TaskRun task, File scriptFile, File envFile, File cmdInputFile, File cmdOutFile ) {
        assert task
        assert scriptFile

        def folder = task.workDirectory

        /*
         * create a script wrapper which do the following
         * 1 - move the TMP directory provided by the sge/oge grid engine
         * 2 - pipe the input stream
         * 3 - launch the user script
         * 4 - un-stage e.g. copy back the result files to the working folder
         */

        def wrapper = new StringBuilder()
        wrapper << 'source ' << envFile.absolutePath << '\n'

        // whenever it has to change to the scratch directory
        def changeDir = changeToScratchDirectory()
        if( changeDir ) {
            wrapper << changeDir << '\n'
        }

        // execute the command script
        wrapper << '('
        if( cmdInputFile ) {
            wrapper << 'cat ' << cmdInputFile << ' | '
        }
        wrapper << "${scriptFile.absolutePath}) &> ${cmdOutFile.absolutePath}" << '\n'

        // "un-stage" the result files
        def resultFiles = taskConfig.outputs.keySet().findAll { it != '-' }
        if( resultFiles && changeDir ) {
            resultFiles.each { name -> wrapper << "for X in $name; do cp \$X $folder; done\n" }
            wrapper << 'rm -rf $NF_SCRATCH &'
        }

        def result = new File(folder, JOB_SCRIPT_FILENAME)
        result.text = task.processor.normalizeScript(wrapper.toString())

        return result
    }

    /**
     * @return The bash script fragment to change to the 'scratch' directory if it has been specified in the task configuration
     */
    protected String changeToScratchDirectory() {

        def scratch = taskConfig.scratch

        if( scratch == null || scratch == false ) {
            return null
        }

        /*
         * when 'scratch' is defined as a bool value
         * try to use the 'TMP' variable, if does not exist fallback to a tmp folder
         */
        if( scratch == true ) {
            return 'NF_SCRATCH=${TMPDIR:-`mktemp -d`} && cd $NF_SCRATCH'
        }

        // convert to string for safety
        scratch = scratch.toString()

        // when it is defined by a variable, just use it
        if( scratch.startsWith('$') ) {
            return "NF_SCRATCH=$scratch && cd \$NF_SCRATCH"
        }

        if( scratch.toLowerCase() in ['ramdisk','ram-disk']) {
            return 'NF_SCRATCH=$(mktemp -d -p /dev/shm/nextflow) && cd $NF_SCRATCH'
        }


        return "NF_SCRATCH=\$(mktemp -d -p $scratch) && cd \$NF_SCRATCH"

    }


    /**
     * Submit the job script to the grid executor
     *
     * @param task The task instance to be executed
     * @param cmdOutFile The file where the job outputs its stdout result
     */
    protected submitJob( TaskRun task, File runnerFile, File cmdOutFile ) {
        assert task

        final folder = task.workDirectory

        // -- log the qsub command
        def cli = getSubmitCommandLine(task)
        log.debug "sub command > '${cli}' -- task: ${task.name}"

        /*
         * launch 'sub' script wrapper
         */
        ProcessBuilder builder = new ProcessBuilder()
                .directory(folder)
                .command( cli as String[] )
                .redirectErrorStream(true)

        // -- configure the job environment
        builder.environment().putAll(task.processor.getProcessEnvironment())

        // -- start the execution and notify the event to the monitor
        Process process = builder.start()
        task.status = TaskRun.Status.RUNNING


        // -- save the 'sub' process output
        def subOutFile = new File(folder, JOB_OUT_FILENAME)
        def subOutStream = new BufferedOutputStream(new FileOutputStream(subOutFile))
        ByteDumper subDumper = new ByteDumper(process.getInputStream(), {  byte[] data, int len -> subOutStream.write(data,0,len) } )
        subDumper.setName("sub_${task.name}")
        subDumper.start()

        try {
            // -- wait the the process completes
            task.exitCode = process.waitFor()
            def success = task.exitCode in taskConfig.validExitCodes
            log.debug "Task completed > ${task.name} -- exit code: ${task.exitCode}; accepted code(s): ${taskConfig.validExitCodes.join(',')}"

            subDumper.await(500)
            subOutStream.close()

            // there may be very loooong delay over NFS, wait at least one minute
            if( success ) {
                Duration.waitFor('60s') { cmdOutFile.exists() }
            }
            if( cmdOutFile.exists() && taskConfig.echo ) {
                print cmdOutFile.text
            }

        }
        finally {
            subDumper.terminate()

            // make sure to release all resources
            IOUtils.closeQuietly(process.in)
            IOUtils.closeQuietly(process.out)
            IOUtils.closeQuietly(process.err)
            process.destroy()

            task.output = collectResultFile(task,'-')
        }

    }

    /**
     * Build up the platform native command line used to submit the job wrapper
     * execution request to the underlying grid, e.g. {@code qsub -q something script.job}
     *
     * @param task The task instance descriptor
     * @return A list holding the command line
     */
    abstract protected List<String> getSubmitCommandLine(TaskRun task)

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
     * Get the file where the task stdout has been saved
     *
     * @param task The task instance for which the output file is required
     * @return The file holding the task stdout
     */
    @Override
    def getStdOutFile( TaskRun task ) {
        assert task

        //  -- return the program output with the following strategy
        //   + program terminated ok -> return the program output output (file)
        //   + program failed and output file not empty -> program output
        //             failed and output EMPTY -> return 'sub' output file
        def cmdOutFile = new File(task.workDirectory, COMMAND_OUTPUT_FILENAME)
        def subOutFile = new File(task.workDirectory, JOB_OUT_FILENAME)
        log.debug "Task cmd output > ${task.name} -- file ${cmdOutFile}; empty: ${cmdOutFile.isEmpty()}"
        log.debug "Task sub output > ${task.name} -- file: ${subOutFile}; empty: ${subOutFile.isEmpty()}"

        def result
        def success = task.exitCode in taskConfig.validExitCodes
        if( success ) {
            result = cmdOutFile.isNotEmpty() ? cmdOutFile : null
        }
        else {
            result = cmdOutFile.isNotEmpty() ? cmdOutFile : ( subOutFile.isNotEmpty() ? subOutFile : null )
        }

        log.debug "Task finished > ${task.name} -- success: ${success}; output: ${result}"
        return result

    }



}
