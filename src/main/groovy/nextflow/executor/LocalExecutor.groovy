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

import java.nio.file.Path

import groovy.util.logging.Slf4j
import nextflow.processor.FileInParam
import nextflow.processor.TaskRun
import nextflow.util.ByteDumper
import nextflow.util.PosixProcess
import org.apache.commons.io.IOUtils
import org.codehaus.groovy.runtime.IOGroovyMethods
/**
 * Executes the specified task on the locally exploiting the underlying Java thread pool
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class LocalExecutor extends AbstractExecutor<LocalProcessHandler> {

    private static final COMMAND_OUT_FILENAME = '.command.out'

    private static final COMMAND_BASH_FILENAME = '.command.sh'

    private static final COMMAND_ENV_FILENAME = '.command.env'

    private static final COMMAND_SCRIPT_FILENAME = '.command.run'


    /**
     * Run a system executable script
     *
     * @param script
     * @return
     */
    @Override
    void launchTask( TaskRun<LocalProcessHandler> task )  {
        assert task
        assert task.@script
        assert task.workDirectory

        final scratch = task.workDirectory
        log.debug "Launching task > ${task.name} -- scratch folder: $scratch"

        /*
         * save the environment to a file
         */
        createEnvironmentFile(task, scratch.resolve(COMMAND_ENV_FILENAME))

        /*
         * save the main script file
         */
        final scriptStr = task.processor.normalizeScript(task.script.toString())
        final scriptFile = scratch.resolve(COMMAND_SCRIPT_FILENAME)
        scriptFile.text = scriptStr
        final interpreter = task.processor.fetchInterpreter(scriptStr)

        /*
         * staging input files
         */
        final files = task.getInputsByType(FileInParam)
        final staging = stagingFilesScript(files)

        /*
         * create the runner script which will launch the script
         */
        def wrapperScript = []
        wrapperScript << '#!/bin/bash -Eeu'
        wrapperScript << "# task: ${task.name}"
        wrapperScript << 'trap onexit 1 2 3 15 ERR'
        wrapperScript << 'function onexit() { local exit_status=${1:-$?}; printf $exit_status > .exitcode; exit $exit_status; }'
        if( staging ) {
            wrapperScript << staging
        }
        wrapperScript << "source $COMMAND_ENV_FILENAME"
        wrapperScript << "$interpreter $COMMAND_SCRIPT_FILENAME"
        wrapperScript << 'onexit'
        wrapperScript << ''

        def wrapperFile = scratch.resolve(COMMAND_BASH_FILENAME)
        wrapperFile.text = task.processor.normalizeScript(wrapperScript.join('\n'))

        // the cmd list to launch it
        List cmd = new ArrayList(taskConfig.shell ?: 'bash' as List ) << COMMAND_BASH_FILENAME
        log.trace "Launch cmd line: ${cmd.join(' ')}"

        /*
         * save the reference to the scriptFile
         */
        task.script = scriptFile

        ProcessBuilder builder = new ProcessBuilder()
                .directory(scratch.toFile())
                .command(cmd)
                .redirectErrorStream(true)

        // the file that will hold the process stdout
        Path outputFile = scratch.resolve(COMMAND_OUT_FILENAME)

        // -- start the execution and notify the event to the monitor
        Process process = builder.start()
        task.stdout = outputFile
        task.handler = new LocalProcessHandler(process, task.name, task.stdin, outputFile, taskConfig.echo)
        task.status = TaskRun.Status.STARTED
        task.submitTimeMillis = task.handler.startTimeMillis
        task.startedTimeMillis = task.handler.startTimeMillis

    }

    @Override
    boolean checkStarted( TaskRun<LocalProcessHandler> task ) {
        return task.isStarted()
    }

    @Override
    boolean checkCompleted( TaskRun<LocalProcessHandler>  task ) {

        if( task.isTerminated() ) {
            return true
        }

        def done = task.handler.hasExited()
        if( done ) {
            task.exitCode = task.handler.exitCode()
            task.status = TaskRun.Status.TERMINATED
            task.handler.destroy()
        }
        else if( taskConfig.maxDuration ) {
            /*
             * check if the task exceed max duration time
             */
            if( task.handler.elapsedTimeMillis() > taskConfig.maxDuration.toMillis() ) {
                task.handler.destroy()
                task.exitCode = task.handler.exitCode()
                task.status = TaskRun.Status.TERMINATED

                // signal has completed
                done = true
            }
        }

        return done
    }



    /**
     * A process wrapper adding the ability to access to the Posix PID
     * and the {@code hasExited} flag
     *
     * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
     */
    static class LocalProcessHandler implements ProcessHandler {

        private final PosixProcess target

        private final startTimeMillis = System.currentTimeMillis()

        private final input

        private final Path output

        private final Boolean echo

        private final String taskName

        private ByteDumper outputDumper

        private BufferedOutputStream outputBuffer

        LocalProcessHandler( Process process, String taskName, def input, Path output, Boolean echo ) {
            this.target = new PosixProcess(process)
            this.taskName = taskName
            this.input = input
            this.output = output
            this.echo = echo

            // handle in/out
            handleInputOutput()
        }

        /**
         * Pipe the process input to the running process and save the stdout to the specified *output* file
         */
        private handleInputOutput() {

            /*
             * Pipe the input data using a parallel thread
             */
            if( input ) {

                Thread.start("Task $taskName > input feeder") {
                    try {
                        IOGroovyMethods.withStream(new BufferedOutputStream(target.getOutputStream())) { writer -> writer << input }
                    }
                    catch( Exception e ) {
                        log.warn "Unable to pipe input data for task: ${taskName}"
                    }
                }
            }

            /*
             * Capture the process stdout and save it to a file
             */
            outputBuffer = output.newOutputStream()
            def writer = { byte[] data, int len ->
                outputBuffer.write(data,0,len)
                if( echo ) { System.out.print(new String(data,0,len)) }
            }

            outputDumper = new ByteDumper(target.getInputStream(), writer)
            outputDumper.setName("Task $taskName > output dumper")
            outputDumper.start()
        }


        long elapsedTimeMillis() {
            System.currentTimeMillis() - startTimeMillis
        }

        /**
         * Check if the submitted job has started
         */
        boolean hasStarted() { target != null }

        /**
         * @return The file containing the job stdout
         */
        Path getOutputFile() { this.output }

        /**
         * Check if the submitted job has terminated its execution
         */
        boolean hasExited() { target.hasExited() }

        /**
         * @return The exit status if the user task
         */
        int exitCode() { target.exitValue() }

        /**
         * Force the submitted job to quit
         */
        void kill() { destroy() }

        /**
         * Destroy the process handler, closing all associated streams
         */
        void destroy() {
            outputDumper.await(500)
            outputDumper.terminate()
            outputBuffer.close()
            IOUtils.closeQuietly(target.getInputStream())
            IOUtils.closeQuietly(target.getOutputStream())
            IOUtils.closeQuietly(target.getErrorStream())
            target.destroy()
        }

    }

}
