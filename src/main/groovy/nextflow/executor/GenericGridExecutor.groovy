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
import groovy.transform.InheritConstructors
import groovy.util.logging.Slf4j
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import nextflow.util.ByteDumper
import nextflow.util.Duration
import org.apache.commons.io.IOUtils

/**
 * Generic task processor executing a task through a grid facility
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@InheritConstructors
abstract class GenericGridExecutor extends ExecutionStrategy {

    protected static final COMMAND_SCRIPT_FILENAME = '.command.sh'

    protected static final COMMAND_OUTPUT_FILENAME = '.command.out'

    protected static final COMMAND_INPUT_FILE = '.command.in'

    protected static final COMMAND_ENV_FILENAME = '.command.env'

    protected static final JOB_OUT_FILENAME = '.job.out'

    protected static final JOB_SCRIPT_FILENAME = '.job.run'


    abstract protected List<String> getSubmitCommandLine(TaskRun task)


    protected String changeToTempFolder() {
        '[ ! -z $TMPDIR ] && cd $TMPDIR'
    }


    @Override
    void launchTask( TaskProcessor processor, TaskRun task ) {
        assert task
        assert task.workDirectory

        final folder = task.workDirectory
        log.debug "Lauching task > ${task.name} -- scratch folder: $folder"

        /*
         * save the environment to a file
         */
        final envMap = processor.getProcessEnvironment()
        final envBuilder = new StringBuilder()
        envMap.each { name, value ->
            if( name ==~ /[a-zA-Z_]+[a-zA-Z0-9_]*/ ) {
                envBuilder << "export $name='$value'" << '\n'
            }
            else {
                log.debug "Task ${task.name} > Invalid environment variable name: '${name}'"
            }
        }
        def envFile = new File(folder, COMMAND_ENV_FILENAME)
        envFile.text = envBuilder.toString()


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
        scriptFile.text = processor.normalizeScript(task.script.toString())
        scriptFile.setExecutable(true)
        task.script = scriptFile

        /*
         * create a script wrapper which do the following
         * 1 - move the TMP directory provided by the sge/oge grid engine
         * 2 - pipe the input stream
         * 3 - launch the user script
         * 4 - un-stage e.g. copy back the result files to the working folder
         */
        File cmdOutFile = new File(folder, COMMAND_OUTPUT_FILENAME)
        def wrapper = new StringBuilder()
        wrapper << 'source ' << envFile.absolutePath << '\n'
        wrapper << changeToTempFolder()  << '\n'

        // execute the command script
        wrapper << '('
        if( cmdInputFile ) {
            wrapper << 'cat ' << cmdInputFile << ' | '
        }
        wrapper << "${scriptFile.absolutePath}) &> ${cmdOutFile.absolutePath}" << '\n'

        // "un-stage" the result files
        def resultFiles = taskConfig.outputs.keySet().findAll { it != '-' }
        if( resultFiles ) {
            wrapper << "if [ \$PWD != $folder ]; then" << '\n'
            resultFiles.each { name -> wrapper << "for X in $name; do cp \$X $folder; done\n" }
            wrapper << 'fi' << '\n'
        }

        new File(folder, JOB_SCRIPT_FILENAME).text = processor.normalizeScript(wrapper.toString())

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
        builder.environment().putAll(processor.getProcessEnvironment())

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
            log.debug "Task completeted > ${task.name} -- exit code: ${task.exitCode}; accepted code(s): ${validExitCodes.join(',')}"

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
