/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package nextflow.executor
import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import nextflow.processor.TaskHandler
import nextflow.processor.TaskId
import nextflow.processor.TaskMonitor
import nextflow.processor.TaskPollingMonitor
import nextflow.processor.TaskRun
import nextflow.scheduler.Protocol.TaskHolder
import nextflow.script.ScriptType
import nextflow.util.Duration
import nextflow.util.ServiceName
/**
 * A Nextflow executor based on Ignite services
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@ServiceName('ignite')
@SupportedScriptTypes( [ScriptType.SCRIPTLET, ScriptType.GROOVY] )
class IgExecutor extends Executor {

    @PackageScope
    IgConnector connector

    /**
     * Initialize the executor by getting a reference to the Ignite connector
     */
    void init() {
        super.init()
        connector = IgConnector.create(taskMonitor)
    }

    /**
     * Creates the task monitor for this executor
     * @return An instance of {@link TaskMonitor}
     */
    @Override
    protected TaskMonitor createTaskMonitor() {
        TaskPollingMonitor.create(session, name, 100, Duration.of('5s'))
    }


    /**
     *  Creates an handler for the specified task
     */
    @Override
    TaskHandler createTaskHandler(TaskRun task) {

        if( task.type == ScriptType.GROOVY ) {
            IgTaskHandler.createGroovyHandler(task, this)
        }
        else {
            IgTaskHandler.createScriptHandler(task, this)
        }

    }


    TaskPollingMonitor getTaskMonitor() {
        (TaskPollingMonitor)super.getTaskMonitor()
    }

    @PackageScope
    void execute( IgBaseTask task ) {
        connector.schedule(task)
    }

    @PackageScope
    boolean checkTaskStarted( TaskId taskId ) {
        connector.checkTaskStarted(taskId)
    }

    @PackageScope
    boolean checkTaskCompleted( TaskId taskId ) {
        connector.checkTaskCompleted(taskId)
    }

    boolean checkTaskFailed( TaskId taskId ) {
        connector.checkTaskFailed(taskId)
    }

    @PackageScope
    TaskHolder removeTaskCompleted( TaskId taskId ) {
        connector.removeTaskCompleted(taskId)
    }

    @PackageScope
    void cancelTask( TaskId taskId ) {
        connector.cancelTask(taskId)
    }

    String dumpQueueStatus() {
        connector.dumpScheduledTasksStatus()
    }

}



