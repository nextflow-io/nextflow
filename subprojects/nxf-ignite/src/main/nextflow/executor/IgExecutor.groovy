/*
 * Copyright (c) 2013-2016, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2016, Paolo Di Tommaso and the respective authors.
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
import nextflow.processor.TaskMonitor
import nextflow.processor.TaskPollingMonitor
import nextflow.processor.TaskRun
import nextflow.script.ScriptType
import nextflow.util.Duration
import org.apache.ignite.IgniteException
import org.apache.ignite.cluster.ClusterNode
import org.apache.ignite.compute.ComputeJob
import org.apache.ignite.compute.ComputeLoadBalancer
import org.apache.ignite.compute.ComputeTaskAdapter
import org.apache.ignite.lang.IgniteCallable
import org.apache.ignite.lang.IgniteFuture
import org.apache.ignite.resources.LoadBalancerResource
import org.jetbrains.annotations.Nullable
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
        TaskPollingMonitor.create(session, name, Duration.of('5s'))
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
    IgniteFuture call( IgniteCallable command ) {
        final compute = connector.compute().withAsync()
        compute.call(command)
        return compute.future()
    }

    @PackageScope
    IgniteFuture execute( ComputeJob task ) {
        final compute = connector.compute().withAsync()
        compute.execute( new IgniteTaskWrapper(task), null)
        compute.future()
    }

    @PackageScope
    boolean checkTaskStarted( taskId ) {
        connector.runningTasks.containsKey(taskId)
    }

    @PackageScope
    void removeTaskCompleted( taskId ) {
        connector.runningTasks.remove(taskId)
    }

    /**
     * An adapter for Ignite compute task
     *
     */
    static class IgniteTaskWrapper extends ComputeTaskAdapter  {

        // Inject load balancer.
        @LoadBalancerResource
        transient ComputeLoadBalancer balancer

        private ComputeJob theJob

        IgniteTaskWrapper( ComputeJob job ) {
            this.theJob = job
        }

        @Override
        Map<? extends ComputeJob, ClusterNode> map(List<ClusterNode> nodes, @Nullable Object arg) throws IgniteException {

            Map<ComputeJob, ClusterNode> jobUnit = [:]
            jobUnit.put(theJob, balancer.getBalancedNode(theJob, null))
            return jobUnit
        }

        @Override
        Object reduce(List list) throws IgniteException {
            return list.get(0)
        }
    }

}



