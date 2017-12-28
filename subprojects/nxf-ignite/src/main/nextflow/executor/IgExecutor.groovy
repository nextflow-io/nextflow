/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
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



