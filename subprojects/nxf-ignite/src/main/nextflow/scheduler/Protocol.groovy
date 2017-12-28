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

package nextflow.scheduler
import groovy.transform.Canonical
import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString
import groovy.util.logging.Slf4j
import nextflow.cloud.CloudDriverFactory
import nextflow.executor.IgBaseTask
import nextflow.processor.TaskId
import nextflow.processor.TaskRun
import nextflow.util.ClusterConfig
import nextflow.util.Duration
import nextflow.util.MemoryUnit
import nextflow.util.SysHelper
import org.apache.ignite.Ignite
/**
 * Defines constants and message objects that implement the Ignite based distributed
 * execution logic
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
interface Protocol {

    static final String PENDING_TASKS_CACHE = 'nextflow.cache.pendingtasks'

    static final String TOPIC_SCHEDULER_EVENTS = 'nextflow.events.scheduler'

    static final String TOPIC_AGENT_EVENTS = 'nextflow.events.agent'

    /**
     * Message sent to notify that a computing node is idle
     */
    @Canonical
    @CompileStatic
    class NodeIdle implements Serializable, Cloneable {

        /**
         * Timestamp when the node entered in the `idle` status
         */
        long idleTimestamp

    }

    /**
     * Message sent from the scheduler to notify the cluster shutdown
     */
    @CompileStatic
    class NodeShutdown implements Cloneable, Serializable {

        public static NodeShutdown INSTANCE = new NodeShutdown()

        private NodeShutdown() {}
    }

    /**
     * Event send from the agent to the scheduler to notify that
     * an instance received a spot/preemptive termination notice
     */
    @Canonical
    @CompileStatic
    class NodeRetired implements Cloneable, Serializable {

        String termination

        String toString() {
            return "termination-notice=$termination"
        }
    }


    /**
     * Message send to all cluster node to notify that a new task has been submitted
     * for execution
     *
     * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
     */
    @CompileStatic
    class TaskAvail {

        /**
         * Singleton object instance
         */
        public static TaskAvail INSTANCE = new TaskAvail()

        private TaskAvail() {}

    }

    /**
     * Message send to a scheduler agent to cancel a job execution
     *
     * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
     */
    @Canonical
    @CompileStatic
    @EqualsAndHashCode
    class TaskCancel implements Serializable, Cloneable {

        /**
         * Id of the task to cancel
         */
        TaskId taskId

        String toString() {
            "taskId=$taskId"
        }
    }
    /**
     * Message send by a scheduler agent to notify that a task execution has completed
     *
     * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
     */
    @CompileStatic
    @EqualsAndHashCode
    class TaskComplete implements Serializable, Cloneable {

        /**
         * The task ID
         */
        TaskId taskId

        /**
         * The result of a task execution either the exit status for a script task
         * or the value returned by a Groovy native task
         */
        def result

        /**
         * The execution thrown during the task execution
         */
        Throwable error

        protected TaskComplete() {}

        static TaskComplete create(IgBaseTask task, result) {
            new TaskComplete(taskId: task.taskId, result: result)
        }

        static TaskComplete error(IgBaseTask task, Throwable error) {
            new TaskComplete(taskId: task.taskId, error: error)
        }

        @Override
        String toString() {
            error ? "taskId=$taskId; error=${error}" : "taskId=$taskId; result=$result"
        }

    }


    /**
     * Holds the metadata of a scheduled task
     *
     * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
     */
    @EqualsAndHashCode
    @CompileStatic
    class TaskHolder implements Serializable, Cloneable {

        /**
         * The actual task object
         */
        IgBaseTask task

        /**
         * The ID of the node where the task is running
         */
        volatile UUID worker

        volatile long submitTimestamp

        volatile long startTimestamp

        volatile boolean started

        volatile boolean completed

        /**
         * Note mark this field as `volatile` because it is accessed by other threads.
         *
         * An alternative implementation would required to define this immutable and
         * have the method {@link #withComplete(nextflow.scheduler.Protocol.TaskComplete)} to create
         * a copy of this object
         *
         * {@link nextflow.executor.IgTaskHandler#checkIfCompleted()}
         */
        volatile Object result

        /**
         * Note mark this field as `volatile` because it is accessed by other threads.
         *
         * An alternative implementation would required to define this immutable and
         * have the method {@link #withComplete(nextflow.scheduler.Protocol.TaskComplete)} to create
         * a copy of this object
         *
         * {@link nextflow.executor.IgExecutor#checkTaskFailed(nextflow.processor.TaskId)}
         */
        volatile Throwable error

        protected TaskHolder() {

        }

        TaskHolder(IgBaseTask task) {
            this.task = task
            submitTimestamp = System.currentTimeMillis()
        }

        TaskHolder withStart( UUID worker ) {
            this.worker = worker
            this.started = true
            this.startTimestamp = System.currentTimeMillis()
            return this
        }

        TaskHolder withComplete( TaskComplete message ) {
            if(!worker) throw new IllegalStateException("Illegal task completion state: `worker` is null -- $this")
            if(!started) throw new IllegalStateException("Illegal task completion state: not started -- $this")

            this.completed = true
            if( message.error )
                this.error = message.error
            else
                this.result = message.result

            return this
        }

        String toString() {
            "taskId=${task.taskId}; worker=$worker; started=$started; completed=$completed; result=$result; error=$error"
        }

        TaskHolder clone() {
            (TaskHolder)super.clone()
        }

        boolean isWaitingMoreThan( Duration duration ) {
            !started && System.currentTimeMillis()-submitTimestamp > duration.toMillis()
        }

    }
    /**
     * Model the computing resources required by a task
     *
     * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
     */
    @ToString
    @EqualsAndHashCode
    @CompileStatic
    class TaskResources implements Serializable, Cloneable {

        int cpus

        MemoryUnit memory

        MemoryUnit disk

        Duration time

        TaskResources() {}

        TaskResources( TaskRun task ) {
            cpus = task.config.getCpus()
            memory = task.config.getMemory()
            disk = task.config.getDisk()
            time = task.config.getTime()
        }

        @Override
        String toString() {
            "cpus=${cpus}; mem=${memory}; disk=${disk}; time=${time}"
        }
    }

    /**
     * Message sent from a scheduler agent to the master node to notify that the task execution has started
     *
     * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
     */
    @Canonical
    @CompileStatic
    class TaskStart implements Serializable, Cloneable {

        /**
         * The ID of the task started
         */
        TaskId taskId

        TaskStart(IgBaseTask task) {
            this.taskId = task.taskId
        }

        protected TaskStart() {}

        @Override
        String toString() {
            "taskId=$taskId"
        }
    }

    @Slf4j
    @CompileStatic
    @EqualsAndHashCode
    class NodeData implements Serializable, Cloneable {

        UUID nodeId

        /**
         * Total amount of resources provided by this node
         */
        Resources resources

        /**
         * Current amount of free resources
         */
        @Deprecated
        Resources free

        /**
         * Host name
         */
        String hostName

        /**
         * Boot time in millis
         */
        long bootTimeMillis

        /**
         * Cloud VM instance ID
         */
        String instanceId

        /**
         * The timestamp when the node become idle
         */
        long idleTimestamp

        static NodeData create(ClusterConfig config) {
            new NodeData(
                    free: new Resources(config),
                    resources: new Resources(config),
                    hostName: SysHelper.getHostName(),
                    bootTimeMillis: SysHelper.getBootTimeMillis() )
        }

        static NodeData create(ClusterConfig config, Ignite ignite) {
            def res = new Resources(config)
            def nodeId = ignite.cluster().localNode().id()
            def instanceId = getCloudInstanceId(config.getCloudDriverName())
            new NodeData(
                    nodeId: nodeId,
                    instanceId: instanceId,
                    free: res,
                    resources: res,
                    hostName: SysHelper.getHostName(),
                    bootTimeMillis: SysHelper.getBootTimeMillis() )
        }


        private static String getCloudInstanceId(String driverName) {
            try {
                return driverName ? CloudDriverFactory.getDriver(driverName).getLocalInstanceId() : null
            }
            catch( Exception e ){
                log.debug "Oops.. Cannot retrieve the instance-id for this node -- May it is not a cloud instance?"
                return null
            }
        }

        Duration getUptime() {
            new Duration(System.currentTimeMillis() - bootTimeMillis)
        }

        boolean isIdle(Duration than) {
            idleTimestamp && System.currentTimeMillis() - idleTimestamp > than.toMillis()
        }

        Duration idle() {
            idleTimestamp ? new Duration(System.currentTimeMillis() - idleTimestamp) : null
        }

        String toString() {
            "nodeId=$nodeId; hostname=$hostName; instance-id: $instanceId; boot-time=${SysHelper.fmtDate(bootTimeMillis)}; tot-res:[$resources]; idle=${idle() ?: '-'}"
        }

    }

    /**
     * Model the resources available for task executions
     */
    @CompileStatic
    @EqualsAndHashCode
    class Resources implements Serializable, Cloneable {

        int cpus

        MemoryUnit memory

        MemoryUnit disk

        protected Resources() {
        }

        Resources(ClusterConfig config) {
            cpus = config.getAttribute('maxCpus', SysHelper.getAvailCpus()) as int
            memory = config.getAttribute('maxMemory', SysHelper.getAvailMemory()) as MemoryUnit
            disk = config.getAttribute('maxDisk', SysHelper.getAvailDisk()) as MemoryUnit
        }

        String toString() {
            "cpus=$cpus; mem=$memory; disk=$disk"
        }


        @Override
        Resources clone() {
            (Resources)super.clone()
        }

    }


}
