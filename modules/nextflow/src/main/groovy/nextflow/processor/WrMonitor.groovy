/*
 * Copyright 2019, Genome Research Limited
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

package nextflow.processor

import java.util.concurrent.LinkedBlockingQueue
import java.util.concurrent.locks.Condition
import java.util.concurrent.locks.Lock
import java.util.concurrent.locks.ReentrantLock
import java.util.concurrent.TimeUnit
import com.google.common.util.concurrent.RateLimiter

import nextflow.processor.TaskFault
import nextflow.processor.BatchContext
import nextflow.processor.BatchHandler
import nextflow.util.Throttle
import nextflow.util.Duration
import nextflow.Session
import nextflow.executor.BatchCleanup
import nextflow.executor.WrRestApi
import nextflow.executor.WrTaskHandler

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j

/**
 * Monitors the queued tasks, adding them to wr and waiting for their termination
 *
 * @author Sendu Bala <sb10@sanger.ac.uk>
 * Based on TaskPollingMonitor by Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * *** currently quite a bit of direct code duplication from TaskPollingMonitor
 */
@Slf4j
@CompileStatic
class WrMonitor implements TaskMonitor {

    private static String RATE_FORMAT = ~/^(\d+\.?\d*)\s*([a-zA-Z]*)/
    final Session session
    final private WrRestApi client
    final Duration dumpInterval

    /**
     * A lock object used to signal the completion of a task execution
     */
    private Lock taskCompleteLock

    /**
     * Locks the queue of pending tasks buffer
     */
    private Lock pendingLock

    /**
     * Condition to signal a new task has been added in the {@link #pendingQueue}
     */
    private Condition taskAvail

    /**
     * A condition that signal when at least a complete task is available
     */
    private Condition taskComplete

    /**
     * Unbounded buffer that holds all tasks scheduled but not yet submitted for execution
     */
    private Queue<WrTaskHandler> pendingQueue

    /**
     * Bounded buffer that holds all {@code TaskHandler}s that have been submitted for execution
     */
    private Queue<TaskHandler> runningQueue

    /**
     * Create the task polling monitor with the provided named parameters object.
     * <p>
     * Valid parameters are:
     * <li>session: The current {@code Session}
     * <li>dumpInterval: Determines how often the executor status is written in the application log file
     *
     * @param params
     */
    protected WrMonitor( Map params ) {
        assert params
        assert params.session instanceof Session
        assert params.client instanceof WrRestApi

        this.session = params.session as Session
        this.client = params.client as WrRestApi
        this.dumpInterval = (params.dumpInterval as Duration) ?: Duration.of('5min')

        this.pendingQueue = new LinkedBlockingQueue()
        this.runningQueue = new LinkedBlockingQueue<TaskHandler>()
    }

    static WrMonitor create( Session session, WrRestApi client ) {
        assert session
        final dumpInterval = session.getMonitorDumpInterval('wr')

        log.debug "Creating task monitor for executor 'wr' > dumpInterval: $dumpInterval"
        new WrMonitor(session: session, client: client, dumpInterval: dumpInterval)
    }

    /**
     * Submits the specified tasks for execution adding them to the queue of scheduled tasks
     *
     * @param handler
     *      A List of {@link WrTaskHandler} instance representing the task to be submitted for execution
     */
    protected void submit(WrTaskHandler handler) {
        // submit the job execution -- throws a ProcessException when submit operation fail
        handler.submit()
        // note: add the 'handler' into the polling queue *after* the submit operation,
        // this guarantees that in the queue are only jobs successfully submitted
        runningQueue.add(handler)
        // notify task submission
        session.notifyTaskSubmit(handler)
    }

    /**
     * Remove a task from task the scheduling queue
     *
     * @param handler
     *      A {@link TaskHandler} instance
     * @return
     *      {@code true} if the task was removed successfully from the tasks polling queue,
     *      {@code false} otherwise
     */
    protected boolean remove(TaskHandler handler) {
        runningQueue.remove(handler)
    }

    /**
     * Schedule a new task for execution
     *
     * @param handler
     *      A {@link TaskHandler} representing the task to be submitted for execution
     */
    @Override
    void schedule(TaskHandler handler) {
        pendingLock.lock()
        try{
            pendingQueue << (WrTaskHandler)handler
            taskAvail.signal()  // signal that a new task is available for execution
            log.trace "Scheduled task > $handler"
        }
        finally {
            pendingLock.unlock()
        }
    }

    /**
     * Evicts a task from the processing tasks queue
     *
     * @param handler
     *      A {@link TaskHandler} instance
     * @return
     *      {@code true} when the specified task was successfully removed from polling queue,
     *      {@code false} otherwise
     */
    @Override
    boolean evict(TaskHandler handler) {
        if( !handler ) {
            return false
        }
        return remove(handler)
    }

    /**
     * Launch the monitoring thread
     *
     * @return
     *      The monitor object itself
     */
    @Override
    TaskMonitor start() {
        log.debug ">>> barrier register (monitor: wr)"
        session.barrier.register(this)

        this.taskCompleteLock = new ReentrantLock()
        this.taskComplete = taskCompleteLock.newCondition()

        this.pendingLock = new ReentrantLock()
        this.taskAvail = pendingLock.newCondition()

        // remove pending tasks on termination
        session.onShutdown { this.cleanup() }

        // launch the thread polling the queue
        Thread.start('Task monitor') {
            try {
                pollLoop()
            }
            finally {
                log.debug "<<< barrier arrives (monitor: wr)"
                session.barrier.arrive(this)
            }
        }

        // launch daemon that submits tasks for execution
        Thread.startDaemon('Task submitter', this.&submitLoop)

        return this
    }

    private void awaitTasks() {
        pendingLock.lock()
        try {
            if( pendingQueue.size()==0 ) {
                taskAvail.await()
            }
        }
        finally {
            pendingLock.unlock()
        }
    }

    /**
     * Wait for new tasks and submit for execution when a slot is available
     */
    protected void submitLoop() {
        while( true ) {
            awaitTasks()
            submitPendingTasks()
        }
    }

    /**
     * Implements the polling strategy
     */
    protected void pollLoop() {

        int iteration=0
        while( true ) {
            final long time = System.currentTimeMillis()
            final tasks = new ArrayList(runningQueue)
            log.trace "Scheduler queue size: ${tasks.size()} (iteration: ${++iteration})"

            // check all running tasks for termination
            checkAllTasks(tasks)

            if( (session.isTerminated() && runningQueue.size()==0 && pendingQueue.size()==0) || session.isAborted() ) {
                break
            }

            await(time)

            if( session.isAborted() ) {
                break
            }

            // dump this line at the dump interval
            Throttle.after(dumpInterval) {
                dumpPendingTasks()
            }
        }
    }

    protected void dumpPendingTasks() {

        try {
            def pending = runningQueue.size()
            if( !pending ) {
                log.debug "No more task to compute -- ${session.dumpNetworkStatus() ?: 'Execution may be stalled'}"
                return
            }

            def msg = []
            msg << "!! executor wr > tasks to be completed: ${runningQueue.size()} -- pending tasks are shown below"
            // dump the first 10 tasks
            def i=0; def itr = runningQueue.iterator()
            while( i++<10 && itr.hasNext() )
                msg << "~> ${itr.next()}"
            if( pending>i )
                msg << ".. remaining tasks omitted."
            log.debug msg.join('\n')
        }
        catch (Throwable e) {
            log.debug "Oops.. expected exception", e
        }
    }

    protected void dumpSubmitQueue() {
        try {
            def pending = pendingQueue.size()
            if( !pending )
                return

            def msg = []
            msg << "%% executor wr > tasks in the submission queue: ${pending} -- tasks to be submitted are shown below"
            // dump the first 10 tasks
            def i=0; def itr = pendingQueue.iterator()
            while( i++<10 && itr.hasNext() )
                msg << "~> ${itr.next()}"
            if( pending>i )
                msg << ".. remaining tasks omitted."
            log.debug msg.join('\n')
        }
        catch (Throwable e) {
            log.debug "Oops.. unexpected exception", e
        }
    }


    /**
     * Await for one or more tasks to be processed
     *
     * @param time
     *      The wait timeout in millis
     */
    protected void await( long time ) {
        def delta = 1000 - (System.currentTimeMillis() - time)
        if( delta <= 0 )
            return

        taskCompleteLock.withLock {
            taskComplete.await( delta, TimeUnit.MILLISECONDS )
        }
    }

    /**
     * Signal that a task has been completed
     */
    @Override
    void signal() {
        taskCompleteLock.withLock {
            taskComplete.signal()
        }
    }

    protected void setupBatchCollector(List<TaskHandler> queue) {
        Map<Class,BatchContext> collectors
        for( int i=0; i<queue.size(); i++ ) {
            final TaskHandler handler = queue.get(i)
            // ignore tasks but BatchHandler
            if( handler instanceof BatchHandler ) {
                // create the main collectors map
                if( collectors == null )
                    collectors = new LinkedHashMap<>()
                // create a collector instance for all task of the same class
                BatchContext c = collectors.getOrCreate(handler.getClass()) { new BatchContext() }
                // set the collector in the handler instance
                handler.batch(c)
            }
        }
    }

    /**
     * Check and update the status of queued tasks
     */
    protected void checkAllTasks(List<TaskHandler> queue) {

        // -- find all task handlers that are *batch* aware
        //    this allows to group multiple calls to a remote system together
        setupBatchCollector(queue)

        // -- iterate over the task and check the status
        for( int i=0; i<queue.size(); i++ ) {
            final handler = queue.get(i)
            try {
                checkTaskStatus(handler)
            }
            catch (Throwable error) {
                handleException(handler, error)
            }
        }

    }

    /**
     * Loop over the queue of pending tasks and submit all
     *
     * @return The number of tasks submitted for execution
     */
    protected int submitPendingTasks() {

        int count = 0
        def itr = pendingQueue.iterator()
        while( itr.hasNext() ) {
            final handler = itr.next()
            try {
                if( session.isSuccess() ) {
                    itr.remove(); count++   // <-- remove the task in all cases
                    submit(handler)
                }
                else
                    break
            }
            catch ( Throwable e ) {
                handleException(handler, e)
                session.notifyTaskComplete(handler)
            }
        }
        log.debug("submitPendingTasks looped through $count pending tasks")
        return count
    }


    final protected void handleException( TaskHandler handler, Throwable error ) {
        def fault = null
        try {
            fault = handler.task.processor.resumeOrDie(handler?.task, error)
        }
        finally {
            // abort the session if a task task was returned
            if (fault instanceof TaskFault) {
                session.fault(fault, handler)
            }
        }
    }


    /**
     * Check the status of the given task
     *
     * @param handler
     *      The {@link TaskHandler} instance of the task to check
     */
    protected void checkTaskStatus( TaskHandler handler ) {
        assert handler

        // check if it is started
        if( handler.checkIfRunning() ) {
            session.notifyTaskStart(handler)
        }

        // check if it is terminated
        if( handler.checkIfCompleted() ) {
            log.debug "Task completed > $handler"

            // since completed *remove* the task from the processing queue
            evict(handler)

            // finalize the tasks execution
            final fault = handler.task.processor.finalizeTask(handler.task)
            // trigger the count down latch when it is a blocking task
            handler.latch?.countDown()

            // notify task completion
            session.notifyTaskComplete(handler)

            // abort the execution in case of task failure
            if (fault instanceof TaskFault) {
                session.fault(fault, handler)
            }
        }

    }

    /**
     * Kill all pending jobs when current execution session is aborted
     */
    protected void cleanup() {
        if( !runningQueue.size() ) return
        log.warn "Killing pending tasks (${runningQueue.size()})"

        def batch = new BatchCleanup()
        while( runningQueue.size() ) {

            TaskHandler handler = runningQueue.poll()
            try {
                // if( handler instanceof WrTaskHandler ) {
                //     ((WrTaskHandler)handler).batch = batch
                // }
                handler.kill()
            }
            catch( Throwable e ) {
                log.debug "Failed to kill pending tasks: ${handler} -- cause: ${e.message}"
            }

            // notify task completion
            handler.task.aborted = true
            session.notifyTaskComplete(handler)
        }

        try {
            batch.kill()
        }
        catch( Throwable e ) {
            log.debug "Failed to kill pending tasks ${batch} -- cause: ${e.message}"
        }
    }
}
