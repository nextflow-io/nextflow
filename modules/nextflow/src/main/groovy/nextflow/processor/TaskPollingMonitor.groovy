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

package nextflow.processor

import java.util.concurrent.ArrayBlockingQueue
import java.util.concurrent.LinkedBlockingQueue
import java.util.concurrent.TimeUnit
import java.util.concurrent.locks.Condition
import java.util.concurrent.locks.Lock
import java.util.concurrent.locks.ReentrantLock

import com.google.common.util.concurrent.RateLimiter
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.executor.BatchCleanup
import nextflow.executor.GridTaskHandler
import nextflow.util.Duration
import nextflow.util.Throttle
/**
 * Monitors the queued tasks waiting for their termination
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@Slf4j
@CompileStatic
class TaskPollingMonitor implements TaskMonitor {

    private static String RATE_FORMAT = ~/^(\d+\.?\d*)\s*([a-zA-Z]*)/
    
    /**
     * The current session object
     */
    final Session session

    /**
     * The time interval (in milliseconds) elapsed which execute a new poll
     */
    final long pollIntervalMillis

    /**
     * Determines how often the executor status is written in the application log file (default: 5min)
     */
    final Duration dumpInterval

    /**
     * The name of the executor associated to this monitor
     */
    final String name


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
     * Condition to signal a new processing slot may be available in the {@link #runningQueue}
     */
    private Condition slotAvail

    /**
     * A condition that signal when at least a complete task is available
     */
    private Condition taskComplete

    /**
     * Unbounded buffer that holds all tasks scheduled but not yet submitted for execution
     */
    private Queue<TaskHandler> pendingQueue

    /**
     * Bounded buffer that holds all {@code TaskHandler}s that have been submitted for execution
     */
    private Queue<TaskHandler> runningQueue

    /**
     * The capacity of the {@code pollingQueue} ie. the max number of tasks can be executed at
     * the same time.
     */
    private int capacity

    /**
     * Define rate limit for task submission
     */
    private RateLimiter submitRateLimit

    /**
     * Create the task polling monitor with the provided named parameters object.
     * <p>
     * Valid parameters are:
     * <li>name: The name of the executor for which the polling monitor is created
     * <li>session: The current {@code Session}
     * <li>capacity: The maximum number of this monitoring queue
     * <li>pollInterval: Determines how often a poll occurs to check for a process termination
     * <li>dumpInterval: Determines how often the executor status is written in the application log file
     *
     * @param params
     */
    protected TaskPollingMonitor( Map params ) {
        assert params
        assert params.session instanceof Session
        assert params.name != null
        assert params.pollInterval != null

        this.name = params.name
        this.session = params.session as Session
        this.pollIntervalMillis = ( params.pollInterval as Duration ).toMillis()
        this.dumpInterval = (params.dumpInterval as Duration) ?: Duration.of('5min')
        this.capacity = (params.capacity ?: 0) as int

        this.pendingQueue = new LinkedBlockingQueue()
        this.runningQueue = capacity ? new ArrayBlockingQueue<TaskHandler>(capacity) : new LinkedBlockingQueue<TaskHandler>()
    }

    static TaskPollingMonitor create( Session session, String name, int defQueueSize, Duration defPollInterval ) {
        assert session
        assert name
        final capacity = session.getQueueSize(name, defQueueSize)
        final pollInterval = session.getPollInterval(name, defPollInterval)
        final dumpInterval = session.getMonitorDumpInterval(name)

        log.debug "Creating task monitor for executor '$name' > capacity: $capacity; pollInterval: $pollInterval; dumpInterval: $dumpInterval "
        new TaskPollingMonitor(name: name, session: session, capacity: capacity, pollInterval: pollInterval, dumpInterval: dumpInterval)
    }

    static TaskPollingMonitor create( Session session, String name, Duration defPollInterval ) {
        assert session
        assert name

        final pollInterval = session.getPollInterval(name, defPollInterval)
        final dumpInterval = session.getMonitorDumpInterval(name)

        log.debug "Creating task monitor for executor '$name' > pollInterval: $pollInterval; dumpInterval: $dumpInterval "
        new TaskPollingMonitor(name: name, session: session, pollInterval: pollInterval, dumpInterval: dumpInterval)
    }

    /**
     * @return The current {@link #runningQueue} instance
     */
    protected Queue<TaskHandler> getRunningQueue() { runningQueue }

    /**
     * @return the current capacity value by the number of slots specified
     */
    int getCapacity() { capacity }


    /**
     * Defines the strategy determining if a task can be submitted for execution.
     *
     * @param handler
     *      A {@link TaskHandler} for a task to be submitted
     * @return
     *      {@code true} if the task satisfies the resource requirements and scheduling strategy implemented
     *      by the polling monitor
     */
    protected boolean canSubmit(TaskHandler handler) {
        (capacity>0 ? runningQueue.size() < capacity : true) && handler.canForkProcess()
    }

    /**
     * Submits the specified task for execution adding it to the queue of scheduled tasks
     *
     * @param handler
     *      A {@link TaskHandler} instance representing the task to be submitted for execution
     */
    protected void submit(TaskHandler handler) {
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
            pendingQueue << handler
            taskAvail.signal()  // signal that a new task is available for execution
            session.notifyTaskPending(handler)
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

        if( remove(handler) ) {
            pendingLock.lock()
            try {
                slotAvail.signal()
                return true
            }
            finally {
                pendingLock.unlock()
            }
        }

        return false
    }

    /**
     * Launch the monitoring thread
     *
     * @return
     *      The monitor object itself
     */
    @Override
    TaskMonitor start() {
        log.trace ">>> barrier register (monitor: ${this.name})"
        session.barrier.register(this)

        this.taskCompleteLock = new ReentrantLock()
        this.taskComplete = taskCompleteLock.newCondition()

        this.pendingLock = new ReentrantLock()
        this.taskAvail = pendingLock.newCondition()
        this.slotAvail = pendingLock.newCondition()

        //
        this.submitRateLimit = createSubmitRateLimit()

        // remove pending tasks on termination
        session.onShutdown { this.cleanup() }

        // launch the thread polling the queue
        Thread.start('Task monitor') {
            try {
                pollLoop()
            }
            finally {
                log.trace "<<< barrier arrives (monitor: ${this.name})"
                session.barrier.arrive(this)
            }
        }

        // launch daemon that submits tasks for execution
        Thread.startDaemon('Task submitter', this.&submitLoop)

        return this
    }

    protected RateLimiter createSubmitRateLimit() {
        def limit = session.getExecConfigProp(name,'submitRateLimit',null) as String
        if( !limit )
            return null

        def tokens = limit.tokenize('/')
        if( tokens.size() == 2 ) {
            /*
             * the rate limit is provide num of task over a duration
             * - eg. 100 / 5 min
             * - ie. max 100 task per 5 minutes
             */

            final X = tokens[0].trim()
            final Y = tokens[1].trim()

            return newRateLimiter(X, Y, limit)
        }

        /*
         * the rate limit is provide as a duration
         * - eg. 200 min
         * - ie. max 200 task per minutes
         */

        final matcher = (limit =~ RATE_FORMAT)
        if( !matcher.matches() )
            throw new IllegalArgumentException("Invalid submit-rate-limit value: $limit -- It must be provide using the following format `num request sec|min|hour` eg. 10 sec ie. max 10 tasks per second")

        final num = matcher.group(1) ?: '_'
        final unit = matcher.group(2) ?: 'sec'

        return newRateLimiter(num, "1 $unit", limit)
    }

    private RateLimiter newRateLimiter( String X, String Y, String limit ) {
        if( !X.isInteger() )
            throw new IllegalArgumentException("Invalid submit-rate-limit value: $limit -- It must be provide using the following format `num request / duration` eg. 10/1s")

        final num = Integer.parseInt(X)
        final duration = Y.isInteger() ? Duration.of( Y+'sec' ) : ( Y[0].isInteger() ? Duration.of(Y) : Duration.of('1'+Y) )
        long seconds = duration.toSeconds()
        if( !seconds )
            throw new IllegalArgumentException("Invalid submit-rate-limit value: $limit -- The interval must be at least 1 second")

        log.debug "Creating submit rate limit of $num reqs by $seconds seconds"
        return RateLimiter.create( num / seconds as double )
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

    private void awaitSlots() {
        pendingLock.lock()
        try {
            slotAvail.await(1, TimeUnit.SECONDS)
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
            // wait for at least at to be available
            awaitTasks()

            // try to submit all pending tasks
            int processed = submitPendingTasks()

            // if no task has been submitted wait for a new slot to be available
            if( !processed ) {
                Throttle.after(dumpInterval) { dumpSubmitQueue() }
                awaitSlots()
            }
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

            // dump this line every two minutes
            Throttle.after(dumpInterval) {
                dumpRunningQueue()
            }
        }
    }

    protected void dumpRunningQueue() {

        try {
            def pending = runningQueue.size()
            if( !pending ) {
                log.debug "No more task to compute -- ${session.dumpNetworkStatus() ?: 'Execution may be stalled'}"
                return
            }

            def msg = []
            msg << "!! executor $name > tasks to be completed: ${runningQueue.size()} -- submitted tasks are shown below"
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
            msg << "%% executor $name > tasks in the submission queue: ${pending} -- tasks to be submitted are shown below"
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
        def delta = this.pollIntervalMillis - (System.currentTimeMillis() - time)
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
     * of which satisfy the {@link #canSubmit(nextflow.processor.TaskHandler)}  condition
     *
     * @return The number of tasks submitted for execution
     */
    protected int submitPendingTasks() {

        int count = 0
        def itr = pendingQueue.iterator()
        while( itr.hasNext() ) {
            final handler = itr.next()
            try {
                submitRateLimit?.acquire()

                if( !canSubmit(handler) )
                    continue

                if( session.isSuccess() ) {
                    itr.remove(); count++   // <-- remove the task in all cases
                    handler.incProcessForks()
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
            // decrement forks count
            handler.decProcessForks()

            // since completed *remove* the task from the processing queue
            evict(handler)

            // finalize the tasks execution
            final fault = handler.task.processor.finalizeTask(handler.task)

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
                if( handler instanceof GridTaskHandler ) {
                    ((GridTaskHandler)handler).batch = batch
                }
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

    /**
     * The queue of tasks pending for submission to the underlying execution system
     */
    protected Queue<TaskHandler> getPendingQueue() {
        return pendingQueue
    }

}

