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
import java.nio.file.Path

import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.processor.DelegateMap
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
@SupportedScriptTypes( [ScriptType.SCRIPTLET, ScriptType.GROOVY] )
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
        TaskPollingMonitor.create(session, name, Duration.of('5s'))
    }


    /**
     *  Creates an handler for the specified task
     */
    @Override
    TaskHandler createTaskHandler(TaskRun task) {

        if( task.type == ScriptType.GROOVY ) {
            return HzTaskHandler.createGroovyHandler(task, taskConfig,this)
        }

        /*
         * otherwise as a bash script
         */
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

        HzTaskHandler.createScriptHandler(task, taskConfig, this)
    }


    TaskPollingMonitor getTaskMonitor() {
        (TaskPollingMonitor)super.getTaskMonitor()
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

    private Path exitFile

    private Path wrapperFile

    private Path outputFile

    private ScriptType type

    /**
     * The member where the job has started
     */
    @PackageScope
    volatile String runningMember

    /**
     * The result object for this task
     */
    @PackageScope
    HzCmdStatus result

    def void setResult( HzCmdStatus result ) {
        this.result = result
        this.lastUpdate = System.currentTimeMillis()
    }

    def void setRunningMember( String value ) {
        this.runningMember = value
        this.lastUpdate = System.currentTimeMillis()
    }


    static HzTaskHandler createScriptHandler( TaskRun task, TaskConfig taskConfig, HzExecutor executor ) {
        def handler = new HzTaskHandler(task,taskConfig)
        handler.executor = executor
        handler.exitFile = task.getCmdExitFile()
        handler.outputFile = task.getCmdOutputFile()
        handler.wrapperFile = task.getCmdWrapperFile()
        handler.type = ScriptType.SCRIPTLET
        return handler
    }

    static createGroovyHandler( TaskRun task, TaskConfig taskConfig, HzExecutor executor ) {
        def handler = new HzTaskHandler(task,taskConfig)
        handler.executor = executor
        handler.type = ScriptType.GROOVY
        return handler
    }

    private HzTaskHandler(TaskRun task, TaskConfig config) {
        super(task,config)
    }

    @Override
    void submit() {

        // clear previous state, eventually
        result = null
        runningMember = null

        // submit to an hazelcast node for execution
        final sender = task.processor.session.uniqueId
        HzCmdCall command
        if( type == ScriptType.SCRIPTLET ) {
            final List cmdLine = new ArrayList(taskConfig.shell ?: 'bash' as List ) << wrapperFile.getName()
            command = new HzCmdCall( sender,  task, cmdLine )
        }
        else {
            command = new HzCmdCall( sender, task )
        }
        executor.connector.executorsQueue.add( command )

        // mark as submitted -- transition to STARTED has to be managed by the scheduler
        status = TaskHandler.Status.SUBMITTED
        log.trace "Task $task > Submitted"
    }

    @Override
    boolean checkIfRunning() {
        if( isSubmitted() && runningMember != null ) {
            log.trace "Task ${task} > RUNNING on member ${runningMember}"
            status = TaskHandler.Status.RUNNING
            return true
        }

        return false
    }

    @Override
    boolean checkIfCompleted() {

        if( isRunning() && result != null && (!exitFile || exitFile.lastModified()>0) ) {
            status = TaskHandler.Status.COMPLETED

            // -- set the task exit code (only when it is a scriptlet task)
            if( isScriptlet() && result.value instanceof Integer )
                task.exitStatus = result.value as Integer

            // -- the task output depend by the kind of the task executed
            if( isScriptlet() )
                task.stdout = outputFile
            else {
                task.stdout = result.value
                task.code.delegate = new DelegateMap( task.processor, result.context )
            }

            // -- set the task result error (if any)
            task.error = result.error

            log.trace "Task ${task} > DONE"
            return true
        }

        return false
    }

    @Override
    void kill() {
        // not implemented
        // see also https://groups.google.com/d/msg/hazelcast/-SVy4k1-QrA/rlUee3i4POAJ
    }


    protected StringBuilder toStringBuilder(StringBuilder builder) {
        super.toStringBuilder(builder)
        builder << "; running: $runningMember"
    }

    /**
     * @return Whenever is a shell script task
     */
    boolean isScriptlet() { type == ScriptType.SCRIPTLET }

    /**
     * @return Whenever is a groovy closure task
     */
    boolean isGroovy() { type == ScriptType.GROOVY }

}

