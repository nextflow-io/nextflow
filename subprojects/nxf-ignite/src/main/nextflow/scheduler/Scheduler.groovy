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
import static nextflow.scheduler.Protocol.PENDING_TASKS_CACHE
import static nextflow.scheduler.Protocol.TOPIC_AGENT_EVENTS
import static nextflow.scheduler.Protocol.TOPIC_SCHEDULER_EVENTS

import javax.cache.CacheException
import java.util.concurrent.BlockingQueue
import java.util.concurrent.ConcurrentHashMap
import java.util.concurrent.LinkedBlockingQueue

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import nextflow.cloud.CloudSpotTerminationException
import nextflow.daemon.IgGridFactory
import nextflow.executor.IgBaseTask
import nextflow.processor.TaskId
import nextflow.processor.TaskPollingMonitor
import nextflow.scheduler.Protocol.NodeData
import nextflow.scheduler.Protocol.NodeIdle
import nextflow.scheduler.Protocol.NodeShutdown
import nextflow.scheduler.Protocol.TaskAvail
import nextflow.scheduler.Protocol.TaskCancel
import nextflow.scheduler.Protocol.TaskComplete
import nextflow.scheduler.Protocol.TaskHolder
import nextflow.scheduler.Protocol.TaskStart
import nextflow.scheduler.Protocol.NodeRetired
import org.apache.ignite.Ignite
import org.apache.ignite.IgniteCache
import org.apache.ignite.IgniteInterruptedException
import org.apache.ignite.events.DiscoveryEvent
import org.apache.ignite.events.Event
import org.apache.ignite.events.EventType
import org.apache.ignite.lang.IgniteBiPredicate
import org.apache.ignite.lang.IgniteCallable
import org.apache.ignite.lang.IgnitePredicate
import org.apache.ignite.resources.IgniteInstanceResource
import org.slf4j.Logger
import org.slf4j.LoggerFactory
/**
 * Implements the scheduler controller logic
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class Scheduler {

    private static final Logger log = LoggerFactory.getLogger(Scheduler)

    static private class ClusterDiscovery implements IgniteCallable<NodeData> {

        @IgniteInstanceResource
        private Ignite ignite

        @Override
        NodeData call() throws Exception {
            def clusterConfig = IgGridFactory.instance().clusterConfig
            return NodeData.create(clusterConfig, ignite)
        }
    }

    /**
     * Local map holding all tasks submitted to a remote agent for execution
     */
    private Map<TaskId,TaskHolder> scheduledTasks

    /**
     * Local map holding completed task until they are not drained by the {@link TaskPollingMonitor}
     */
    private Map<TaskId,TaskHolder> completedTasks

    /**
     * Distributed map holding all task to be executed
     */
    private IgniteCache<TaskId,IgBaseTask> pendingTasks

    /**
     * Map each node UUID to the associated host name
     */
    private Map<UUID, NodeData> workerNodes

    /**
     * {@link Ignite} instance
     */
    private Ignite ignite

    /**
     * Reference to the {@link TaskPollingMonitor} instance
     */
    private TaskPollingMonitor monitor

    /**
     * Holds a queue of received message to be processed
     */
    private BlockingQueue<Closure> messageQueue = new LinkedBlockingQueue<>()

    /**
     * Thread which process the messages in the {@link #messageQueue}
     */
    private Thread messageProcessor

    /**
     * The auto-scaling policy
     */
    private Autoscaler autoscaler

    private boolean cloudEnabled

    private long bootTimestamp


    /**
     * Initialize the scheduler instance
     *
     * @param ignite An {@link Ignite} instance
     * @param monitor A {@link TaskPollingMonitor} instance
     */
    Scheduler init(Ignite ignite, TaskPollingMonitor monitor) {
        assert ignite
        assert monitor

        this.ignite = ignite
        this.monitor = monitor
        this.bootTimestamp = System.currentTimeMillis()
        this.workerNodes = new ConcurrentHashMap<>()
        this.scheduledTasks = new ConcurrentHashMap<>()
        this.completedTasks = new ConcurrentHashMap<>()
        this.pendingTasks = ignite.cache(PENDING_TASKS_CACHE)

        discoverWorkers()
        createEventProcessor()
        registerEvents()

        return this
    }

    void registerAutoscaler( Autoscaler autoscaler ) {
        this.autoscaler = autoscaler
        this.autoscaler.init(workerNodes, scheduledTasks)
        this.cloudEnabled = true
    }

    /**
     * Discover the nodes that made-up the cluster
     */
    private void discoverWorkers() {

        def nodes = ignite.compute().broadcast(new ClusterDiscovery())

        def buffer = new StringBuilder("+++ Initial cluster topology:\n")
        nodes.each { node ->
            workerNodes[node.nodeId] = node
            buffer << '- ' << node.toString() << '\n'
        }

        log.debug buffer.toString()
    }


    @PackageScope boolean isRunning() {
        messageProcessor.isAlive()
    }

    private void createEventProcessor() {

        messageProcessor = Thread.start('scheduler-thread') {

            while( true ) {
                try {
                    messageQueue.take().call()
                }
                catch( InterruptedException e ) {
                    // time to die
                    break
                }
                catch( Throwable e ) {
                    log.debug("+++ Can't process received message", e)
                }
            }
        }
    }

    /**
     * Dispatch a generic message to a concrete message handler
     *
     * @param sender The {@link UUID} identifier of the node sending the message
     * @param message The message object e.g. {@link TaskStart}
     */
    private IgniteBiPredicate<UUID,Object> createMessageDispatcher() {

        { UUID sender, Object message ->

            //
            // note: the message is not processed inline, but added as a closure
            // in the `messageQueue` object to be processed orderly in the same thread
            //
            messageQueue << {

                if( message instanceof TaskStart ) {
                    onTaskStart(sender, message)
                }
                else if( message instanceof TaskComplete ) {
                    onTaskComplete(sender, message)
                }
                else if( message instanceof NodeData ) {
                    onNodeStart(sender, message)
                }
                else if( message instanceof NodeIdle ) {
                    onNodeIdle(sender, message)
                }
                else if( message instanceof NodeRetired ) {
                    onNodeRetired(sender, message)
                }
                else {
                    throw new IllegalArgumentException("Unknown worker message: $message")
                }

            }
            return true

        } as IgniteBiPredicate<UUID,Object>

    }

    /**
     * Dispatch a generic event to a concrete message handler
     *
     * @param event The {@link Event} object received
     */
    private IgnitePredicate<Event> createEventDispatcher() {

        { Event event ->

            messageQueue << {
                if( event instanceof DiscoveryEvent ) {
                    if( event.type() == EventType.EVT_NODE_LEFT ) {
                        onNodeLeft(event.eventNode().id())
                        return
                    }
                    if( event.type() == EventType.EVT_NODE_FAILED ) {
                        onNodeFailed(event.eventNode().id())
                        return
                    }

                }

                throw new IllegalArgumentException("Unknown event: $event")
            }
            return true

        } as IgnitePredicate<Event>

    }

    /**
     * Registers the events and messages listen by the scheduler
     *
     * Note: events received are appended to the {@link #messageQueue} to be
     * processed asynchronously. This is required otherwise Ignite may stall (hang)
     * when sending a new message in the same thread of the arriving message
     */
    private registerEvents() {
        // -- listen for events from the scheduler agent to this class
        ignite
            .message()
            .localListen( TOPIC_SCHEDULER_EVENTS, createMessageDispatcher() )

        // -- listen for a node left event
        ignite
            .events()
            .localListen( createEventDispatcher(), EventType.EVT_NODE_LEFT )

        // -- listen for a node left failed
        ignite
            .events()
            .localListen( createEventDispatcher(), EventType.EVT_NODE_FAILED )
    }

    /**
     * Schedule one or more tasks for execution. Each {@link IgBaseTask} object is
     * added to the {@link #pendingTasks} distributed cache from where it will be
     * picked from a remote scheduler agent to be processed
     *
     * @param tasks One or more {@link IgBaseTask} instances
     *
     */
    void schedule( IgBaseTask... tasks ) {
        messageQueue << { schedule0(tasks) }
    }

    private void schedule0( IgBaseTask... tasks ) {

        log.trace "+++ Scheduling tasks: taskId=${tasks.collect{ IgBaseTask t -> t.taskId }.join(',')}"

        for( int i=0; i<tasks.size(); i++ ) {
            final task = tasks[0]
            // before add table of scheduled tasks
            scheduledTasks.put(task.getTaskId(), new TaskHolder(task))
            // after append to the queue of pending tasks -- this will trigger the execution on remote workers
            pendingTasks.put(task.getTaskId(), task)
        }

        // notify workers a new task is available
        notifyTaskAvail()
    }

    /**
     * Send a message {@link TaskAvail} to the remote nodes
     */
    private void notifyTaskAvail() {
        ignite.message().send( TOPIC_AGENT_EVENTS, TaskAvail.INSTANCE )
    }

    private String hostName(UUID nodeId) {
        def node = workerNodes.get(nodeId)
        if ( !node ) return null
        cloudEnabled ? node.instanceId : node.hostName
    }

    private void onNodeStart(UUID sender, NodeData message) {
        assert sender == message.nodeId
        workerNodes[sender] = message
        autoscaler?.onNodeStart(message)
        // note: log *after* updating the `workerNodes` map otherwise `hostName` will report a null
        log.debug "+++ Node joined the cluster: [${hostName(sender)}] $sender"
    }

    private void onNodeIdle(UUID sender, NodeIdle message ) {
        def node = workerNodes[sender]
        if( !node ) {
            log.debug "+++ Unknown idle node: nodeId=$sender"
            return
        }

        if( message.idleTimestamp < this.bootTimestamp ) {
            // this means that the work was idle before the computation started
            // just ignore it
            return
        }

        node.idleTimestamp = message.idleTimestamp
        log.debug "+++ Node idle: [${hostName(sender)}] not working for ${node.idle()}"
    }


    /**
     * Handler invoked when a {@link TaskStart} is received
     *
     * @param sender The ID of the node that sent the message
     * @param message A {@link TaskStart} object which identifies the task started
     */
    private void onTaskStart(UUID sender, TaskStart message) {

        def holder = scheduledTasks.get(message.taskId)
        if( holder ) {
            log.trace "+++ Task started: $message [${hostName(sender)}] $sender"
            holder.withStart(sender)
            // -- reset the idle attribute, if any
            def node = workerNodes.get(sender)
            if( node ) {
                node.idleTimestamp = 0
            }
            else {
                log.debug "+++ Oops.. Can't find node data with id=$sender"
            }
        }
        else {
            log.debug "+++ Oops.. Started task is unknown -- [${hostName(sender)}] $message"
        }
    }

    /**
     * Handler invoked when a {@link TaskComplete} message is received
     * @param sender The ID of the node that sent the message
     * @param message A {@link TaskComplete} object which identifies the task complete along with its result
     */
    private void onTaskComplete(UUID sender, TaskComplete message) {

        def holder = scheduledTasks.get(message.taskId)
        if( holder ) {
            log.trace "+++ Task complete: $message [${hostName(sender)}] $sender"
            completedTasks.put(message.taskId, holder.withComplete(message))
            scheduledTasks.remove(message.taskId)
            // -- notify that a task has completed
            monitor.signal()
        }
        else {
            log.debug "+++ Got task complete message but cannot find it -- $message [${hostName(sender)}]"
        }
    }

    /**
     * Handle invoked when a node leave the cluster. The tasks that where assigned
     * to that node are rescheduled to be executed on a different node.
     *
     * @param nodeId The id of the cluster node that has left the cluster
     */
    private void onNodeFailed(UUID nodeId) {
        log.debug "+++ Node failed: [${hostName(nodeId)}] $nodeId"
        removeRunningTaskOnNode(nodeId,'failed')
        workerNodes.remove(nodeId)
    }

    /**
     * Method handler invoked when a cluster node fail is detected. The tasks
     * assigned to the failing node are rescheduled to be executed on a different node.
     *
     * @param nodeId The id of the cluster node that fail
     */
    private void onNodeLeft(UUID nodeId) {
        log.debug "+++ Node left: [${hostName(nodeId)}] $nodeId"
        removeRunningTaskOnNode(nodeId,'leaving')
        workerNodes.remove(nodeId)
    }

    /**
     * Method handler invoked when a cluster node received a spot/preemptive termination
     * noticed e.g. it is going to be retired by the provide and terminated in a few seconds
     *
     * @param nodeId
     */
    private void onNodeRetired(UUID nodeId, NodeRetired message) {
        log.debug "+++ Node retired: [${hostName(nodeId)}] $nodeId -- $message"
        removeRunningTaskOnNode(nodeId, 'retired')
        workerNodes.remove(nodeId)
    }

    private void removeRunningTaskOnNode(UUID nodeId, String reason) {
        // -- finds all tasks whose worker id matches the specified `nodeId`
        def tasks = (List<IgBaseTask>)scheduledTasks
                .values()
                .findResults { TaskHolder it -> it.worker == nodeId ? it.task : null }

        // -- reschedule matching tasks for execution
        if( !tasks ) {
            log.trace "+++ No pending task on $reason node: [${hostName(nodeId)}]"
            return
        }

        log.trace "+++ Dropping tasks on $reason node: [${hostName(nodeId)}] taskId=${tasks.collect{ it.taskId }.join(', ') ?: 'n/a'}"
        def itr = tasks.iterator()
        while( itr.hasNext() ) {
            def task = itr.next()
            // -- simulate an error message
            def cause = (reason=='retired'
                        ? new CloudSpotTerminationException("Computing node was retired: [${hostName(nodeId)}]")
                        : new RuntimeException("Task aborted due to failure on node: [${hostName(nodeId)}]") )
            def failure = TaskComplete.error(task, cause)
            onTaskComplete(nodeId, failure)
        }
    }

    /**
     * Check if the task with the specified id has started
     *
     * @param taskId The task identifier
     * @return {@code true} if the task has started, or {@code false} otherwise
     */
    boolean checkTaskStarted( TaskId taskId ) {
        scheduledTasks.get(taskId)?.started || completedTasks.containsKey(taskId)
    }

    /**
     * Check if the task with the specified id has completed
     *
     * @param taskId The task identifier
     * @return {@code true} if the task has completed, or {@code false} otherwise
     */
    boolean checkTaskCompleted( TaskId taskId ) {
        completedTasks.containsKey(taskId)
    }

    boolean checkTaskFailed( TaskId taskId ) {
        completedTasks.get(taskId)?.error != null
    }

    /**
     * Cancel the execution of a task
     *
     * @param taskId The {@link TaskId} ID of the task execution to cancel
     */
    void cancelTask( TaskId taskId ) {

        messageQueue << {

            log.trace "+++ Cancelling task: taskId=${taskId}"
            boolean removed = false
            try {
                 removed = pendingTasks.remove(taskId)
            }
            catch (CacheException e) {
                if( !(e.cause instanceof IgniteInterruptedException) )
                    throw e
            }

            def holder = scheduledTasks.get(taskId)
            if( holder ) {
                if( holder.worker ) {
                    def worker = ignite.cluster().forNodeId(holder.worker)
                    ignite.message(worker).send( TOPIC_AGENT_EVENTS, new TaskCancel(taskId) )
                }
                scheduledTasks.remove(taskId)
            }

            if( !removed && !holder ) {
                log.trace "+++ Oops.. Unable to cancel task: taskId=${taskId}"
            }
        }
    }

    /**
     * Get the task runtime information and remove from the {@link #scheduledTasks} structure
     *
     * @param taskId The id of the task to retrieve
     * @return A {@link TaskHolder} containing the task runtime info or {@code null} if the task is not available
     */
    TaskHolder removeTaskCompleted( TaskId taskId ) {
        def result = completedTasks.get(taskId)
        completedTasks.remove(taskId)
        return result
    }

    String dumpScheduledTasksStatus() {
        def result = new StringBuilder()
        def itr = scheduledTasks.values().iterator()
        while( itr.hasNext() ) {
            result << itr.next().toString()
        }
        return result.toString()
    }

    /**
     * Shutdown scheduler remote agents by sending a {@link Protocol#TOPIC_AGENT_EVENTS} message
     */
    void shutdownRemoteAgents() {
        def nodes = ignite.cluster().forRemotes()
        ignite.message(nodes).send( TOPIC_AGENT_EVENTS, NodeShutdown.INSTANCE )
    }

    /**
     * Shutdown the scheduler object
     */
    void shutdownScheduler() {
        messageProcessor?.interrupt()
        autoscaler?.closeQuietly()
    }
}
