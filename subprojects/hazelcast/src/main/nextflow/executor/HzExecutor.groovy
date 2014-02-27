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

import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.processor.TaskConfig
import nextflow.processor.TaskHandler
import nextflow.processor.TaskMonitor
import nextflow.processor.TaskPollingMonitor
import nextflow.processor.TaskRun
import nextflow.script.ScriptType
import nextflow.util.Duration
/**
 * Creates an executor that submits process to an Hazelcast cluster for execution
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class HzExecutor extends AbstractExecutor  {

    @PackageScope
    HzConnector connector

    /**
     * Initialize the executor by getting a reference to the Hazelcast connector
     */
    def void init() {
        super.init()
        connector = HzConnector.create(taskMonitor)
    }

    /**
     * Creates the task monitor for this executor
     * @return An instance of {@code TaskMonitor}
     */
    @Override
    protected TaskMonitor createTaskMonitor() {
        TaskPollingMonitor.create(session, name, Duration.of('1s'))
    }


    /**
     *  Creates an handler for the specified task
     */
    @Override
    TaskHandler createTaskHandler(TaskRun task) {

        if( task.type == ScriptType.GROOVY ) {
            throw new UnsupportedOperationException("Native tasks are not supported")
        }

        /*
         * otherwise as a bash script
         */
        final bash = new BashWrapperBuilder(task)
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

        new HzTaskHandler(task, taskConfig, this)
    }


    TaskPollingMonitor getTaskMonitor() {
        (TaskPollingMonitor)super.getTaskMonitor()
    }


    /**
     * The sender ID is the current session ID
     */
    @PackageScope
    UUID getSenderId() {
        return session.uniqueId
    }

}



/**
 * A task handler for Hazelcast cluster
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class HzTaskHandler extends TaskHandler {

    private HzExecutor executor

    private final Path exitFile

    private final Path wrapperFile

    private final Path outputFile

    @PackageScope
    volatile HzCmdResult result

    protected HzTaskHandler(TaskRun task, TaskConfig taskConfig, HzExecutor executor) {
        super(task, taskConfig)
        this.executor = executor
        this.exitFile = task.getCmdExitFile()
        this.outputFile = task.getCmdOutputFile()
        this.wrapperFile = task.getCmdWrapperFile()
    }

    @Override
    void submit() {
        final List cmd = new ArrayList(taskConfig.shell ?: 'bash' as List ) << wrapperFile.getName()

        // submit to an hazelcast node for execution
        def command = new HzCmdCall( executor.senderId, task.id, task.workDirectory, cmd )
        executor.connector.executorsQueue.add( command )

        // mark as submitted -- transition to STARTED has to be managed by the scheduler
        status = TaskHandler.Status.SUBMITTED
    }

    @Override
    boolean checkIfRunning() {
        if( isSubmitted() ) {
            log.trace "Task ${task.name} > RUNNING"
            status = TaskHandler.Status.RUNNING
            return true
        }

        return false
    }

    @Override
    boolean checkIfCompleted() {

        if( isRunning() && result != null ) {
            // TODO Add a check on the result file existence due to NFS latency
            status = TaskHandler.Status.COMPLETED

            // -- set the task exit code (only when it is a scriptlet task)
            if( result.isScriptlet() && result.value instanceof Integer )
                task.exitStatus = result.value as Integer

            // -- the task output depend by the kind of the task executed
            if( result.isScriptlet() )
                task.stdout = outputFile
            else {
                task.stdout = result.value
                task.code.delegate = result.context
            }

            // -- set the task result error (if any)
            task.error = result.error

            log.trace "Task ${task.name} > DONE (${task.exitStatus})"
            return true
        }

        return false
    }

    @Override
    void kill() {
        // not implemented
        // see also https://groups.google.com/d/msg/hazelcast/-SVy4k1-QrA/rlUee3i4POAJ
    }

}

