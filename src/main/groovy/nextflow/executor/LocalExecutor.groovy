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
import nextflow.processor.TaskConfig
import nextflow.processor.TaskHandler
import nextflow.processor.TaskRun
import nextflow.util.PosixProcess
import org.apache.commons.io.IOUtils
import org.codehaus.groovy.runtime.IOGroovyMethods

/**
 * Executes the specified task on the locally exploiting the underlying Java thread pool
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class LocalExecutor extends AbstractExecutor {

    @Override
    def LocalTaskHandler createTaskHandler(TaskRun task) {

        assert task
        assert task.@script
        assert task.workDirectory

        log.debug "Launching task > ${task.name} -- work folder: ${task.workDirectory}"

        final bash = new BashWrapperBuilder(task)

        // set the environment
        bash.environment = task.processor.getProcessEnvironment()
        bash.environment.putAll( task.getInputEnvironment() )

        // staging input files
        bash.stagingScript = {
            final files = task.getInputsByType(FileInParam)
            final staging = stagingFilesScript(files)
            return staging
        }

        // create the wrapper script
        bash.build()

        new LocalTaskHandler( task, taskConfig )
    }


}


/**
 * A process wrapper adding the ability to access to the Posix PID
 * and the {@code hasExited} flag
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class LocalTaskHandler extends TaskHandler {

    private final startTimeMillis = System.currentTimeMillis()

    private final Path exitFile

    private final Long maxDurationMillis

    private final Path wrapperFile

    private final Path outputFile

    private PosixProcess process

    private boolean destroyed

    LocalTaskHandler( TaskRun task, TaskConfig taskConfig ) {
        super(task, taskConfig)
        // create the task handler
        this.exitFile = task.getCmdExitFile()
        this.outputFile = task.getCmdOutputFile()
        this.wrapperFile = task.getCmdWrapperFile()
        this.maxDurationMillis = taskConfig.maxDuration?.toMillis()
    }

    @Override
    def void submit() {

        // the cmd list to launch it
        List cmd = new ArrayList(taskConfig.shell ?: 'bash' as List ) << wrapperFile.getName()
        log.trace "Launch cmd line: ${cmd.join(' ')}"

        ProcessBuilder builder = new ProcessBuilder()
                .directory(task.workDirectory.toFile())
                .command(cmd)
                .redirectErrorStream(true)

        // -- start the execution and notify the event to the monitor
        process = new PosixProcess(builder.start())

        // handle in/out
        pipeStdInput()

        // mark as STARTED
        status = Status.RUNNING
    }

    /**
     * Pipe the process input to the running process and save the stdout to the specified *output* file
     */
    private pipeStdInput() {

        /*
         * Pipe the input data using a parallel thread
         */
        final input = task.stdin
        if( !input ) { return }

        Thread.start("Task ${task.name} > input feeder") {
            try {
                IOGroovyMethods.withStream(new BufferedOutputStream(process.getOutputStream())) { writer -> writer << input }
            }
            catch( Exception e ) {
                log.warn "Unable to pipe input data for task: ${task.name}"
            }
        }

    }


    long elapsedTimeMillis() {
        System.currentTimeMillis() - startTimeMillis
    }

    /**
     * Check if the submitted job has started
     */
    @Override
    boolean checkIfRunning() {
        process != null && status == Status.RUNNING
    }

    /**
     * Check if the submitted job has terminated its execution
     */
    @Override
    boolean checkIfTerminated() {

        if( status == Status.TERMINATED ) {
            return true
        }

        def done = process.hasExited()
        if( done ) {
            task.exitCode = process.exitValue()
            task.stdout = outputFile
            status = Status.TERMINATED
            destroy()
            return true
        }

        if( maxDurationMillis ) {
            /*
             * check if the task exceed max duration time
             */
            if( elapsedTimeMillis() > maxDurationMillis ) {
                destroy()
                task.exitCode = process.exitValue()
                task.stdout = outputFile
                status = Status.TERMINATED

                // signal has completed
                return true
            }
        }

        return false
    }


    /**
     * Force the submitted job to quit
     */
    @Override
    void kill() { destroy() }

    /**
     * Destroy the process handler, closing all associated streams
     */
    void destroy() {

        if( destroyed ) { return }

        IOUtils.closeQuietly(process.getInputStream())
        IOUtils.closeQuietly(process.getOutputStream())
        IOUtils.closeQuietly(process.getErrorStream())
        process.destroy()
        destroyed = true
    }

}
