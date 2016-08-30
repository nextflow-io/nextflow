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

package nextflow.scheduler
import static nextflow.Const.ROLE_MASTER
import static nextflow.scheduler.Protocol.PENDING_TASKS_CACHE
import static nextflow.scheduler.Protocol.TOPIC_AGENT_EVENTS
import static nextflow.scheduler.Protocol.TOPIC_SCHEDULER_EVENTS

import java.nio.channels.ClosedByInterruptException
import java.util.concurrent.BlockingQueue
import java.util.concurrent.Callable
import java.util.concurrent.ConcurrentHashMap
import java.util.concurrent.ExecutorService
import java.util.concurrent.Executors
import java.util.concurrent.Future
import java.util.concurrent.LinkedBlockingQueue
import java.util.concurrent.RejectedExecutionException
import java.util.concurrent.TimeUnit
import java.util.concurrent.locks.Condition
import java.util.concurrent.locks.Lock
import java.util.concurrent.locks.ReentrantLock

import groovy.transform.CompileDynamic
import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.transform.TupleConstructor
import groovy.util.logging.Slf4j
import nextflow.cloud.CloudDriver
import nextflow.cloud.CloudDriverFactory
import nextflow.cloud.CloudSpotTerminationException
import nextflow.daemon.IgGridFactory
import nextflow.executor.IgBaseTask
import nextflow.processor.TaskId
import nextflow.scheduler.Protocol.NodeData
import nextflow.scheduler.Protocol.NodeIdle
import nextflow.scheduler.Protocol.NodeShutdown
import nextflow.scheduler.Protocol.Resources
import nextflow.scheduler.Protocol.TaskAvail
import nextflow.scheduler.Protocol.TaskCancel
import nextflow.scheduler.Protocol.TaskComplete
import nextflow.scheduler.Protocol.TaskStart
import nextflow.util.ClusterConfig
import nextflow.util.Duration
import nextflow.util.MemoryUnit
import nextflow.util.SysHelper
import org.apache.ignite.Ignite
import org.apache.ignite.IgniteCache
import org.apache.ignite.cache.query.ScanQuery
import org.apache.ignite.cluster.ClusterGroupEmptyException
import org.apache.ignite.events.DiscoveryEvent
import org.apache.ignite.events.Event
import org.apache.ignite.events.EventType
import org.apache.ignite.lang.IgniteBiPredicate
import org.apache.ignite.lang.IgnitePredicate
/**
 * Implements the scheduler execution logic. Each worker deploy an instance of this class
 * to process the tasks submitted by driver nextflow application
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class SchedulerAgent implements Closeable {


    /**
     * Predicate that scans all tasks whose resource request match the ones available in the scheduler agent
     */
    static class MatchingResources implements IgniteBiPredicate<TaskId,IgBaseTask> {

        int cpus
        MemoryUnit memory
        MemoryUnit disk

        MatchingResources( Resources avail ) {
            cpus = avail.cpus
            memory = avail.memory
            disk = avail.disk
        }

        @Override
        boolean apply(TaskId taskId, IgBaseTask task) {
            if(task.resources.cpus > cpus) return false
            if(task.resources.memory && task.resources.memory > memory) return false
            if(task.resources.disk && task.resources.disk > disk) return false

            return true
        }

        String toString() {
            "cpus=$cpus; mem=$memory; disk=$disk"
        }
    }

    @CompileStatic
    private class Processor extends Thread {

        private Lock checkpoint = new ReentrantLock()

        private Condition notEmpty = checkpoint.newCondition()

        private Resources current

        private volatile boolean stopped

        private BlockingQueue<Closure> messageQueue = new LinkedBlockingQueue<>()

        private long idleTimestamp

        private long _1_min = Duration.of('1 min').toMillis()

        Processor() {
            this.name = 'scheduler-agent'
        }

        @Override
        void run() {
            this.current = new Resources(config)
            log.debug "=== Agent resources: $current"

            while( !stopped ) {
                try {
                    if( masterId ) {
                        processMessages()
                        checkSpotTermination()
                        if( processPendingTasks(current) ) continue
                        checkIfIdle()
                        // wait for new messages
                        checkpoint.withLock {
                            notEmpty.await(5, TimeUnit.SECONDS)
                        }
                    }
                    else if( !stopped ) {
                        resetState()
                        waitForMasterNodeToJoin()
                    }

                }
                catch( InterruptedException e ) {
                    log.debug("=== Message processor interrupted")
                    stopped = true
                }
                catch( RejectedExecutionException e ) {
                    log.debug "=== Task execution rejected -- ${e.message ?: e}"
                }
                catch( CloudSpotTerminationException e ) {
                    log.debug "=== Detected spot termination notice: $e.termination -- Starting shutdown"
                    stopped = true
                    abortPendingTasks(e)
                    close(true)
                }
                catch( Exception e ) {
                    log.error "=== Unexpected scheduler agent error", e
                }
            }

            // consume latest messages
            while( messageQueue ) {
                processMessages()
            }
        }

        private void abortPendingTasks(Exception e) {
            if( runningTasks ) {
                log.debug "=== aborting pending tasks: taskId=${runningTasks.keySet().join(",") ?: 'n/a'}"

                def itr = runningTasks.values().iterator()
                while( itr.hasNext() ) {
                    RunHolder holder = itr.next()
                    holder.exception = e
                    holder.future.cancel(true)
                }
            }

        }

        private void checkSpotTermination() {
            def termination = driver?.getLocalTerminationNotice()
            if( termination )
                throw new CloudSpotTerminationException(termination)

            if( simulateSpotTermination && runningTasks )
                throw new CloudSpotTerminationException('fake-termination')
        }

        private void resetState() {

            if( runningTasks ) {
                log.debug "=== Cancelling running tasks: taskId=${runningTasks.keySet().join(', ') ?: 'n/a'}"
                def itr = runningTasks.values().iterator()
                while( itr.hasNext() ) {
                    RunHolder holder = itr.next()
                    holder.future.cancel(true)
                }
                runningTasks.clear()
            }

            // reset the distributed cache
            pendingTasks.clear()

            current = new Resources(config)
            log.debug "=== Agent resources after reset: $current"
            messageQueue.clear()
            idleTimestamp = 0
        }

        private void waitForMasterNodeToJoin() {
            int c=0
            while( !masterId ) {
                if ( c++ % 60 == 0 ) {
                    log.debug "=== Waiting for master node to join.."
                }

                try {
                    sleep 5_000
                }
                catch (InterruptedException e) {
                    stopped = true
                }
            }
        }

        void shutdown() {
            stopped = true
            newMessage()
        }

        void newMessage() {
            checkpoint.withLock {
                notEmpty.signal()
            }
        }

        void sendMessage( Closure closure ) {
            messageQueue << closure
            newMessage()
        }

        /**
         * Process the messages added to the {@link #messageQueue} queue
         *
         * @param res
         */
        void processMessages() {
            def Closure msg
            while( (msg=messageQueue.poll()) ) {
                msg.call()
            }
        }


        int processPendingTasks( Resources avail ) {
            def tasks = pendingTasks
                    .query(new ScanQuery<TaskId,IgBaseTask>(new MatchingResources(avail)))
                    .getAll()
                    .collect { it -> it.value }

            def count=0
            def itr = tasks.iterator()
            while( itr.hasNext() && avail.cpus && avail.memory ) {
                count++
                final it = itr.next()
                final res = it.resources

                if( !canRun(it, avail) )
                    continue

                if( pendingTasks.getAndRemove(it.taskId) ) {
                    log.debug "=== Picked task up: taskId=${it.taskId}"
                    // -- decrement resources
                    avail.cpus -= res.cpus
                    avail.memory -= res.memory
                    avail.disk -= res.disk

                    // -- send to the executor
                    try {
                        runningTasks[it.taskId] = new RunHolder(executor.submit( runTask(it) ))
                    }
                    catch (RejectedExecutionException e) {
                        rollbackResources(it)
                        throw e
                    }

                    // -- reset the idle status
                    idleTimestamp = 0

                    // -- signal that task has started
                    notifyTaskStart(it)
                }
            }

            return count
        }

        boolean canRun( IgBaseTask it, Resources avail ) {

            final res = it.resources

            if( res.cpus && res.cpus > avail.cpus ) {
                log.debug "=== Cannot execute task: taskId=${it.taskId} -- CPUs request exceed available (req=${res.cpus}; avail=${avail.cpus})"
                return false
            }

            if( res.memory && res.memory > avail.memory ) {
                log.debug "=== Cannot execute task: taskId=${it.taskId} -- Memory request exceed available (req=${res.memory}; avail=${avail.memory})"
                return false
            }

            if( res.disk && res.disk > avail.disk ) {
                log.debug "=== Cannot execute task: taskId=${it.taskId} -- Disk request exceed available (req=${res.disk}; avail=${avail.disk})"
                return false
            }

            return true
        }

        Callable runTask(IgBaseTask task) {
            return { runTask0(task) } as Callable
        }

        void runTask0( IgBaseTask task ) {

            try {
                def result = task.call()
                notifyComplete(task, result)
            }
            catch( InterruptedException | ClosedByInterruptException e ) {
                def holder = runningTasks[task.taskId]
                if( holder?.exception ) {
                    log.debug "=== Task execution was interrupted: taskId=${task.taskId} -- Cause: $holder.exception"
                    notifyErrorAsync(task, holder.exception)
                }
                else {
                    log.debug "=== Task execution was interrupted: taskId=${task.taskId} -- Message: ${e.message ?: e}"
                }
            }
            catch( Throwable error ) {
                notifyError(task, error)
            }
            finally {
                // Note: since this method is run concurrently by an executor service
                // invoke the `restoreResources` method by using the message to the owning thread
                // In this way it not necessary to apply an expensive synchronisation
                // logic when increasing/decreasing resources amount
                sendMessage({ owner.rollbackResources(task) })
            }

        }

        @CompileDynamic
        void notifyErrorAsync( IgBaseTask task, Exception e) {
            Thread.start { notifyError(task,e) }
        }

        /**
         * Restore the resources "consumed" by the task
         *
         * @param task The {@link IgBaseTask} that used some resources to be restored
         */
        void rollbackResources(IgBaseTask task) {
            final used = task.resources
            final taskId = task.taskId

            if( runningTasks.containsKey(taskId)) {
                runningTasks.remove(taskId)
                // restore resources
                current.cpus += used.cpus
                current.memory += used.memory
                current.disk = SysHelper.getAvailDisk()
                log.debug "=== Resources after task execution: taskId=$taskId; $current"

                // track the time when enter in idle state
                if( current.cpus == total.cpus ) {
                    idleTimestamp = System.currentTimeMillis()
                }
            }
        }

        void checkIfIdle() {

            final now = System.currentTimeMillis()
            if( idleTimestamp && now-idleTimestamp >_1_min ) {
                // send a message to notify the `idle` status of this node
                notifyNodeIdle(idleTimestamp)
                // reset the timestamp to avoid to send multiple notifications for the same idle condition
                idleTimestamp = 0
            }
        }
    }

    /**
     * Pair object holding a {@link Future} to a running task
     */
    @TupleConstructor
    private static class RunHolder {
        Future future
        volatile Exception exception
    }

    /**
     * The underlying executor service
     */
    private ExecutorService executor

    /**
     * Distributed cached of all tasks waiting to be processed
     */
    private IgniteCache<TaskId,IgBaseTask> pendingTasks

    /**
     * Local map of the current running tasks
     */
    private Map<TaskId,RunHolder> runningTasks = new ConcurrentHashMap<>()

    /**
     * Reference to the {@link Ignite} instance
     */
    private Ignite ignite

    /**
     * Reference to the {@link ClusterConfig} object
     */
    private ClusterConfig config

    private Processor processor

    private Resources total

    private volatile boolean closed

    private volatile UUID masterId

    private CloudDriver driver

    private boolean simulateSpotTermination

    /**
     * Initialise the scheduler agent
     *
     * @param ignite The {@link Ignite} instance
     * @param config The {@link ClusterConfig} instance
     */
    SchedulerAgent(Ignite ignite, ClusterConfig config, UUID masterId = null) {
        this.config = config
        this.ignite = ignite
        this.pendingTasks = ignite.cache(PENDING_TASKS_CACHE)
        this.executor = Executors.newFixedThreadPool(SysHelper.getAvailCpus())
        this.total = new Resources(config)
        this.driver = getCloudDriver(config)
        this.processor = new Processor()
        this.simulateSpotTermination = config.getAttribute('simulateSpotTermination') as boolean

        // -- register events to listen to
        registerEvents()

        // -- notify the node has started
        this.masterId = masterId ?: getMasterNodeId()
        if( this.masterId ) {
            notifyNodeStart()
        }
    }

    private CloudDriver getCloudDriver( ClusterConfig config ) {
        final driverName = config.getCloudDriverName()
        try {
            return config.isCloudCluster() ? CloudDriverFactory.get(driverName) : null
        }
        catch( Exception e ) {
            log.error "=== Can't load cloud driver: `$driverName`", e
            return null
        }
    }

    private IgniteBiPredicate<UUID,Object> createMessageDispatcher() {

        //
        // agent messages
        //
        return { UUID uuid, Object message ->

            if( message instanceof TaskAvail ) {
                processor.newMessage()
            }
            else if( message  instanceof TaskCancel ) {
                onCancelTask(message)
            }
            else if( message instanceof NodeShutdown ) {
                onNodeShutdown(uuid)
            }
            else {
                throw new IllegalStateException("Unknown agent event: ${message?.getClass()?.getName()}")
            }
            return true

        } as IgniteBiPredicate<UUID,Object>

    }

    private IgnitePredicate createEventDispatcher() {

        return { Event event ->

            if( event instanceof DiscoveryEvent ) {
                if( event.type() == EventType.EVT_NODE_LEFT ) {
                    onNodeLeft(event.eventNode().id())
                }
                else if( event.type() == EventType.EVT_NODE_FAILED ) {
                    onNodeFailed(event.eventNode().id())
                }
                else if( event.type() == EventType.EVT_NODE_JOINED ) {
                    onNodeJoined(event.eventNode().id())
                }
                else
                    throw new IllegalArgumentException("Unknown event: $event")
            }

            return true

        } as IgnitePredicate<Event>

    }

    private void onNodeJoined(UUID nodeId) {
        if( !masterId && isMasterNode(nodeId) ) {
            log.debug "=== Master node joined: nodeId=$nodeId"
            masterId = nodeId
        }
        else {
            log.debug "=== Cluster node joined: nodeId=$nodeId"
        }
    }

    private void onNodeFailed(UUID nodeId) {
        if( nodeId == masterId ) {
            log.debug "=== Master node failed: nodeId=$nodeId"
            masterId = null
            processor.newMessage()
        }
    }

    private void onNodeLeft(UUID nodeId) {
        if( nodeId == masterId ) {
            log.debug "=== Master node left: nodeId=$nodeId"
            masterId = null
            processor.newMessage()
        }
    }

    private void onNodeShutdown(UUID nodeId) {
        close(true)
    }


    /**
     * Register event the scheduler agent listen for
     */
    private void registerEvents() {

        ignite.message().localListen(TOPIC_AGENT_EVENTS, createMessageDispatcher())

        def dispatcher = createEventDispatcher()
        ignite
                .events()
                .localListen( dispatcher, EventType.EVT_NODE_FAILED )

        ignite
                .events()
                .localListen( dispatcher, EventType.EVT_NODE_LEFT )

        ignite
                .events()
                .localListen( dispatcher, EventType.EVT_NODE_JOINED )

    }

    /**
     * Launch the scheduler agent execution logic. It checks for task in
     * the {@link #pendingTasks} distributed cached, pick the ones which resources
     * match the avail ones and execute them.
     *
     * @return The {@link SchedulerAgent} instance itself
     */
    SchedulerAgent run() {
        processor.start()
        return this
    }

    private isMasterNode(UUID nodeId) {
        ignite.cluster().node(nodeId).attribute(IgGridFactory.NODE_ROLE) == ROLE_MASTER
    }

    private UUID getMasterNodeId() {
        ignite.cluster().forAttribute(IgGridFactory.NODE_ROLE, ROLE_MASTER)?.node()?.id()
    }

    private void sendMessageToMaster( String topic, message ) {
        if( !masterId ) {
            log.debug "=== Master node is unknown -- Cannot send message: [${message.getClass().getSimpleName()}] $message"
            return
        }

        try {
            final master = ignite .cluster() .forNodeId(masterId)
            ignite .message(master) .sendOrdered(topic, message, 0)
        }
        catch ( ClusterGroupEmptyException e ) {
            log.debug "=== Master node is not available -- Cannot send message: [${message.getClass().getSimpleName()}] $message"
        }
    }

    /**
     * Notify the scheduler that the task execution has started by sending a
     * {@link TaskStart} message
     */
    @PackageScope void notifyTaskStart(IgBaseTask task) {
        sendMessageToMaster(TOPIC_SCHEDULER_EVENTS, new TaskStart(task))
    }

    @PackageScope void notifyNodeStart() {
        def data = NodeData.create(config, ignite)
        sendMessageToMaster(TOPIC_SCHEDULER_EVENTS, data)
    }

    /**
     * Notify the scheduler that a task execution has completed by sending a {@link TaskComplete} message
     *
     * @param task
     * @param result
     */
    @PackageScope void notifyComplete(IgBaseTask task, result) {
        try {
            log.debug "=== Notify task complete: taskId=${task.taskId}; result=$result"
            final payload = TaskComplete.create(task, result)
            sendMessageToMaster(TOPIC_SCHEDULER_EVENTS, payload)
        }
        catch( Exception e ) {
            log.error "=== Failed to notify task completion: taskId=${task.taskId}; result=$result", e
        }
    }

    /**
     * Notify the scheduler that a task execution has failed by sending a {@link TaskComplete} message
     *
     * @param task
     * @param error
     */
    @PackageScope void notifyError(IgBaseTask task, Throwable error) {
        try {
            final taskId = task.taskId
            log.debug "=== Notify task complete [error]: taskId=${taskId}; error=$error"
            final payload = TaskComplete.error(task, error)
            sendMessageToMaster(TOPIC_SCHEDULER_EVENTS, payload)
        }
        catch( Exception e ) {
            log.error "=== Failed to notify task completion: taskId=${task.taskId}; error=$error", e
        }
    }

    @PackageScope void notifyNodeIdle(long last) {
        log.debug "=== Notify node idle"
        sendMessageToMaster(TOPIC_SCHEDULER_EVENTS, new NodeIdle(last))
    }

    /**
     * Method handler invoked when receiving a {@link TaskCancel} message. It
     * cancels the specified running tasks.
     *
     * @param message
     */
    @PackageScope void onCancelTask( TaskCancel message ) {
        def holder = runningTasks.get(message.taskId)
        if( holder ) {
            log.debug "=== Cancelling task: taskId=${message.taskId}"
            holder.future.cancel(true)
        }
        else {
            log.debug "=== Unable to find task to cancel: taskId=${message.taskId}"
        }
    }

    /**
     * Shutdown the scheduler agent
     */
    synchronized void close(boolean shutdownIgnite = false) {
        log.debug "=== Scheduler agent shutting down"
        if( closed ) return

        // -- shutdown executor
        executor.shutdownNow()
        processor.shutdown()

        // -- shutdown the ignite instance
        if( !shutdownIgnite )
            return

        Thread.start {
            print "Cleaning up .. "
            sleep 3_000;    // give a few seconds to send pending messages
            ignite.close()
            println "Done."
        }

        closed = true
    }


}
