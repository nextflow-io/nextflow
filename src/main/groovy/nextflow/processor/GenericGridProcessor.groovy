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

package nextflow.processor
import groovy.transform.InheritConstructors
import groovy.util.logging.Slf4j
import nextflow.util.ByteDumper
import nextflow.util.Duration
import nextflow.util.MemoryUnit
import org.apache.commons.io.IOUtils
/**
 * Generic task processor executing a task throgh a grid facility
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@InheritConstructors
abstract class GenericGridProcessor extends AbstractTaskProcessor {

    protected static final COMMAND_SCRIPT_FILENAME = '.command.sh'

    protected static final COMMAND_OUTPUT_FILENAME = '.command.out'

    protected static final COMMAND_INPUT_FILE = '.command.input'

    protected static final JOB_OUT_FILENAME = '.job.out'

    protected static final JOB_SCRIPT_FILENAME = '.job.sh'

    protected String queue

    protected MemoryUnit maxMemory

    protected Duration maxDuration

    /**
     * The SGE/OGE cluster queue to which submit the job
     */

    GenericGridProcessor queue( String queue0 ) {
        this.queue = queue0
        return this
    }

    /**
     * The max memory allow to be used to the job, this value will set the 'virtual_free' qsub cli option
     * <p>
     * Read more about SGE virtual_free vs mem_free at the following links
     * http://gridengine.org/pipermail/users/2011-December/002215.html
     * http://www.gridengine.info/tag/virtual_free/
     *
     * @param mem0 The maximum amount of memory expressed as string value,
     *              accepted units are 'B', 'K', 'M', 'G', 'T', 'P'. So for example
     *              {@code maxMemory '100M'}, {@code maxMemory '2G'}, etc.
     */
    GenericGridProcessor maxMemory( String mem0 ) {
        this.maxMemory = new MemoryUnit(mem0)
        return this
    }

    /**
     * The max duration time allowed for the job to be executed, this value sets the '-l h_rt' squb command line option.
     *
     * @param duration0 The max allowed time expressed as duration string, Accepted units are 'min', 'hour', 'day'.
     *                  For example {@code maxDuration '30 min'}, {@code maxDuration '10 hour'}, {@code maxDuration '2 day'}
     */
    GenericGridProcessor maxDuration( String duration0 ) {
        this.maxDuration = new Duration(duration0)
        return this
    }


    abstract protected List<String> getSubmitCommandLine(TaskRun task)


    @Override
    protected void launchTask(TaskRun task) {
        assert task
        assert task.workDirectory

        final folder = task.workDirectory
        log.debug "Lauching task > ${task.name} -- scratch folder: $folder"

        /*
         * save the command input (if any)
         * the file content will be piped to the executed user script
         */
        File cmdInputFile = null
        if( task.input != null ) {
            cmdInputFile = new File(folder, COMMAND_INPUT_FILE)
            cmdInputFile.text = task.input
        }

        /*
         * save the 'user' script to be executed
         */
        def scriptFile = new File(folder, COMMAND_SCRIPT_FILENAME)
        scriptFile.text = task.script.toString()
        // keep a reference to the file
        task.script = scriptFile

        /*
         * create a script wrapper which do the following
         * 1 - move the TMP directory provided by the sge/oge grid engine
         * 2 - pipe the input stream
         * 3 - launch the user script
         * 4 - un-stage e.g. copy back the result files to the working folder
         */
        File cmdOutFile = new File(folder, COMMAND_OUTPUT_FILENAME)
        def scriptShell = shell instanceof List ? shell.join(' ') : shell
        def wrapper = new StringBuilder()
        wrapper << '#!/bin/bash -ue' << '\n'
        wrapper << '[ ! -z $TMP ] && cd $TMP' << '\n'
        wrapper << '[ ! -z $TMPDIR ] && cd $TMPDIR'  << '\n'

        // execute the command script
        wrapper << '('
        if( cmdInputFile ) {
            wrapper << 'cat ' << cmdInputFile << ' | '
        }
        wrapper << "$scriptShell $scriptFile) &> ${cmdOutFile.absolutePath}" << '\n'

        // "un-stage" the result files
        def resultFiles = outputs.keySet().findAll { it != '-' }
        if( resultFiles ) {
            wrapper << "if [ \$PWD != $folder ]; then" << '\n'
            resultFiles.each { name -> wrapper << "for X in $name; do cp \$X $folder; done\n" }
            wrapper << 'fi' << '\n'
        }

        new File(folder, JOB_SCRIPT_FILENAME).text = wrapper.toString()

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
        builder.environment().putAll(getProcessEnvironment())

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
            def success = task.exitCode in validExitCodes
            log.debug "Task completeted > ${task.name} -- exit code: ${task.exitCode}; accepted code(s): ${validExitCodes.join(',')}"

            subDumper.await(500)
            subOutStream.close()

            // there may be very loooong delay over NFS, wait at least one minute
            if( success ) {
                Duration.waitFor('60s') { cmdOutFile.exists() }
            }
            if( cmdOutFile.exists() && echo ) {
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



    protected getStdOutFile( TaskRun task ) {
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
        def success = task.exitCode in validExitCodes
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
