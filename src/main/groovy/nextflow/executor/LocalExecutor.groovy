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
import org.apache.commons.io.IOUtils
import org.codehaus.groovy.runtime.IOGroovyMethods

/**
 * Executes the specified task on the locally exploiting the underlying Java thread pool
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class LocalExecutor extends AbstractExecutor {

    private static final COMMAND_OUT_FILENAME = '.command.out'

    private static final COMMAND_RUNNER_FILENAME = '.command.run'

    private static final COMMAND_ENV_FILENAME = '.command.env'

    private static final COMMAND_SCRIPT_FILENAME = '.command.sh'


    /**
     * Run a system executable script
     *
     * @param script
     * @return
     */
    @Override
    void launchTask( TaskRun task )  {
        assert task
        assert task.@script
        assert task.workDirectory

        final scratch = task.workDirectory
        log.debug "Lauching task > ${task.name} -- scratch folder: $scratch"

        /*
         * save the environment to a file
         */
        saveEnvironment(task, new File(scratch,COMMAND_ENV_FILENAME))

        /*
         * save the main script file
         */
        def scriptFile = new File(scratch, COMMAND_SCRIPT_FILENAME)
        scriptFile.text = task.processor.normalizeScript(task.script.toString())
        scriptFile.setExecutable(true)

        /*
         * create the runner script which will launch the script
         */
        def runnerText = """
                    source ${COMMAND_ENV_FILENAME}
                    ./${COMMAND_SCRIPT_FILENAME}
                    """
        def runnerFile = new File(scratch, COMMAND_RUNNER_FILENAME)
        runnerFile.text = task.processor.normalizeScript(runnerText)
        runnerFile.setExecutable(true)

        /*
         * save the reference to the scriptFile
         */
        task.script = scriptFile

        ProcessBuilder builder = new ProcessBuilder()
                .directory(scratch)
                .command(runnerFile.absolutePath)
                .redirectErrorStream(true)

        // -- start the execution and notify the event to the monitor
        Process process = builder.start()
        task.status = TaskRun.Status.RUNNING

        // -- copy the input value to the process standard input
        if( task.input != null ) {
            pipeTaskInput( task, process )
        }

        File fileOut = new File(scratch, COMMAND_OUT_FILENAME)
        ByteDumper dumper = null
        try {
            // -- print the process out if it is not capture by the output
            //    * The byte dumper uses a separate thread to capture the process stdout
            //    * The process stdout is captured in two condition:
            //      when the flag 'echo' is set or when it goes in the output channel (outputs['-'])
            //
            BufferedOutputStream streamOut = new BufferedOutputStream( new FileOutputStream(fileOut) )

            def handler = { byte[] data, int len ->
                streamOut.write(data,0,len)
                if( taskConfig.echo ) System.out.print(new String(data,0,len))
            }
            dumper = new ByteDumper(process.getInputStream(), handler)
            dumper.setName("dumper-$taskConfig.name")
            dumper.start()

            // -- wait the the process completes
            if( taskConfig.maxDuration ) {
                log.debug "Running task > ${task.name} -- waiting max: ${taskConfig.maxDuration}"
                process.waitForOrKill(taskConfig.maxDuration.toMillis())
                task.exitCode = process.exitValue()
            }
            else {
                log.debug "Running task > ${task.name} -- wait forever"
                task.exitCode = process.waitFor()
            }

            log.debug "Task completed > ${task.name} -- exit code: ${task.exitCode}; success: ${task.exitCode in taskConfig.validExitCodes}"

            dumper?.await(500)
            streamOut.close()

        }
        finally {

            dumper?.terminate()
            IOUtils.closeQuietly(process.in)
            IOUtils.closeQuietly(process.out)
            IOUtils.closeQuietly(process.err)
            process.destroy()

            task.output = fileOut
        }
    }

    @Override
    def getStdOutFile( TaskRun task ) {

        new File(task.workDirectory, COMMAND_OUT_FILENAME)

    }


    /**
     * Pipe the {@code TaskDef#input} to the {@code Process}
     *
     * @param task The current task to be executed
     * @param process The system process that will run the task
     */
    protected void pipeTaskInput( TaskRun task, Process process ) {

        Thread.start {
            try {
                IOGroovyMethods.withStream(new BufferedOutputStream(process.getOutputStream())) {writer -> writer << task.input}
            }
            catch( Exception e ) {
                log.warn "Unable to pipe input data for task: ${task.name}"
            }
        }
    }

}
