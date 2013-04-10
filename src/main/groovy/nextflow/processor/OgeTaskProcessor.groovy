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
 * Execute a task script by running it on the SGE/OGE cluster
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@InheritConstructors
class OgeTaskProcessor extends AbstractTaskProcessor {

    private static final COMMAND_SCRIPT_FILENAME = '.command.sh'

    private static final COMMAND_OUTPUT_FILENAME = '.command.out'

    private static final QSUB_SCRIPT_FILENAME = '.qsub.sh'

    private static final QSUB_OUT_FILENAME = '.qsub.out'

    private String queue

    private String qsubCmdLine

    private MemoryUnit maxMemory

    private Duration maxDuration

    /**
     * The SGE/OGE cluster queue to which submit the job
     */

    OgeTaskProcessor queue( String queue0 ) {
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
    OgeTaskProcessor maxMemory( String mem0 ) {
        this.maxMemory = new MemoryUnit(mem0)
        return this
    }

    /**
     * The max duration time allowed for the job to be executed, this value sets the '-l h_rt' squb command line option.
     *
     * @param duration0 The max allowed time expressed as duration string, Accepted units are 'min', 'hour', 'day'.
     *                  For example {@code maxDuration '30 min'}, {@code maxDuration '10 hour'}, {@code maxDuration '2 day'}
     */
    OgeTaskProcessor maxDuration( String duration0 ) {
        this.maxDuration = new Duration(duration0)
        return this
    }

    /**
     * Extra options appended to the generated 'qsub' command line
     */
    OgeTaskProcessor qsubCmdLine( String cmdLine ) {
        this.qsubCmdLine = cmdLine
        return this
    }

    /*
     * Prepare the 'qsub' cmdline. The following options are used
     * - wd: define the job working directory
     * - terse: output just the job id on the output stream
     * - o: define the file to which redirect the standard output
     * - e: define the file to which redirect the error output
     */
    protected List<String> getQsubCommandLine(TaskDef task) {

        final result = []

        result << 'qsub'
        result << '-terse'
        result << '-wd $PWD'
        result << '-o' << COMMAND_OUTPUT_FILENAME
        result << '-j y'
        result << '-sync y'
        result << '-V'


        // add other parameters (if any)
        if(queue) {
            result << '-q'  << queue
        }

        if( maxDuration ) {
            result << "-l h_rt=${maxDuration.format('HH:mm:ss')}"
        }

        if( maxMemory ) {
            result << "-l virtual_free=${maxMemory.toString().replaceAll(/[\sB]/,'')}"
        }

        // -- the job name
        result << '-N' << "nf-${name}-${task.index}"

        // -- at the end append the command script wrapped file name
        if ( qsubCmdLine ) {
            result << qsubCmdLine
        }

        // -- last entry to 'script' file name
        result << COMMAND_SCRIPT_FILENAME

        return result
    }


    @Override
    protected void launchTask(TaskDef task) {
        assert task
        assert task.workDirectory

        final scratch = task.workDirectory
        log.debug "Lauching task > ${task.name} -- scratch folder: $scratch"

        /*
         * save the original script to be executed
         */
        def scriptFile = new File(scratch, COMMAND_SCRIPT_FILENAME).absoluteFile
        scriptFile.text = task.script.toString()

        // -- keep a reference to the file
        task.script = scriptFile

        // -- log the qsub command
        def cli = getQsubCommandLine(task).join(' ')
        log.debug "qsub command > '${cli}' -- task: ${task.name}"

        /*
         * Save the 'qsub' command line
         */
        new File(scratch, QSUB_SCRIPT_FILENAME).text = cli

        /*
         * launch 'qsub' script wrapper
         */
        ProcessBuilder builder = new ProcessBuilder()
                .directory(scratch)
                .command( shell, QSUB_SCRIPT_FILENAME )
                .redirectErrorStream(true)

        // -- configure the job environment
        builder.environment().putAll(getProcessEnvironment())

        // -- start the execution and notify the event to the monitor
        Process process = builder.start()
        task.status = TaskDef.Status.RUNNING

        // -- pipe the input value to the process standard input
        if( task.input != null ) {
            try {
                process.withOutputStream{ writer -> writer << task.input }
            }
            catch( IOException e ) {
                log.warn "Unable to pipe input data for task: ${task.name}", e
            }
        }

        // -- save the 'qsub' process output
        def qsubOutFile = new File(scratch, QSUB_OUT_FILENAME)
        def qsubOutStream = new BufferedOutputStream(new FileOutputStream(qsubOutFile))
        ByteDumper qsubDumper = new ByteDumper(process.getInputStream(), {  byte[] data, int len -> qsubOutStream.write(data,0,len) } )
        qsubDumper.setName("qsub-$name")
        qsubDumper.start()

        // -- print the process out if it is not capture by the output
        //    * The byte dumper uses a separate thread to capture the process stdout
        //    * The process stdout is captured in two condition:
        //      when the flag 'echo' is set or when it goes in the output channel (outputs['-'])
        //
        File cmdOutFile = new File(scratch, COMMAND_OUTPUT_FILENAME)
        def handler = echo ? { byte[] data, int len ->  System.out.print(new String(data,0,len)) } : null
        ByteDumper cmdDumper = new ByteDumper( cmdOutFile, handler )
        cmdDumper.setName("dumper-$name")
        cmdDumper.start()

        try {
            // -- wait the the process completes
            task.exitCode = process.waitFor()
            def success = task.exitCode in validExitCodes
            log.debug "Task completeted > ${task.name} -- exit code: ${task.exitCode}; accepted code(s): ${validExitCodes.join(',')}"

            qsubDumper.await(500)
            qsubOutStream.close()

            // there may be very loooong delay over NFS, wait at least one minute
            if(success) {
                cmdDumper.await(60_000)
            }

            // save the exit-code
            new File(scratch, '.exitcode').text = task.exitCode

        }
        finally {
            qsubDumper.terminate()
            cmdDumper.terminate()

            // make sure to release all resources
            IOUtils.closeQuietly(process.in)
            IOUtils.closeQuietly(process.out)
            IOUtils.closeQuietly(process.err)
            process.destroy()

            task.output = collectResultFile(task,'-')
        }

    }


    protected collectResultFile(TaskDef task, String fileName ) {
        assert fileName
        assert task

        if( fileName == '-' ) {

            //  -- return the program output with the following strategy
            //   + program terminated ok -> return the program output output (file)
            //   + program failed and output file not empty -> program output
            //             failed and output EMPTY -> return 'qsub' output file
            def cmdOutFile = new File(task.workDirectory, COMMAND_OUTPUT_FILENAME)
            def qsubOutFile = new File(task.workDirectory, QSUB_OUT_FILENAME)
            log.debug "Task cmd output > ${task.name} -- file ${cmdOutFile}; empty: ${cmdOutFile.isEmpty()}"
            log.debug "Task qsub output > ${task.name} -- file: ${qsubOutFile}; empty: ${qsubOutFile.isEmpty()}"

            def result
            def success = task.exitCode in validExitCodes
            if( success ) {
                result = cmdOutFile.isNotEmpty() ? cmdOutFile : null
            }
            else {
                result = cmdOutFile.isNotEmpty() ? cmdOutFile : ( qsubOutFile.isNotEmpty() ? qsubOutFile : null )
            }

            log.debug "Task finished > ${task.name} -- success: ${success}; output: ${result}"
            return result

        }
        else {
            super.collectResultFile(task,fileName)
        }
    }



}
