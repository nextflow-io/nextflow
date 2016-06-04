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

package nextflow.processor
import java.util.concurrent.ConcurrentLinkedQueue
import java.util.concurrent.TimeUnit
import java.util.concurrent.locks.Condition
import java.util.concurrent.locks.Lock
import java.util.concurrent.locks.ReentrantLock

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.executor.BatchCleanup
import nextflow.executor.GridTaskHandler
import nextflow.util.Duration
/**
 *
 * Monitors the queued tasks waiting for their termination
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@Slf4j
@CompileStatic
class TaskPollingMonitor implements TaskMonitor {

    /**
     * The current session object
     */
    final Session session

    /**
     * The tasks dispatcher
     */
    final TaskDispatcher dispatcher

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
     * A lock object used to access in a synchronous manner the polling queue
     */
    private Lock tasksQueueLock

    /**
     * A lock object used to signal the completion of a task execution
     */
    private Lock taskCompleteLock

    /**
     * A condition that signal when at least a complete task is available
     */
    private Condition taskComplete

    private Condition notFull

    private Queue<TaskHandler> pollingQueue

    private int capacity

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
        this.dispatcher = session.dispatcher
        this.pollIntervalMillis = ( params.pollInterval as Duration ).toMillis()
        this.dumpInterval = (params.dumpInterval as Duration) ?: Duration.of('5min')
        this.capacity = (params.capacity ?: 0) as int

        this.pollingQueue = new ConcurrentLinkedQueue<>()
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

    protected Queue<TaskHandler> getPollingQueue() { pollingQueue }

    public TaskDispatcher getDispatcher() { dispatcher }

    /**
     * Set the number of tasks the monitor will handle in parallel.
     * This value must be changed by using a mutex lock
     *
     * @param value
     */
    final void capacitySet( int newValue ) {
        if( newValue > capacity ) {
            capacityInc(newValue-capacity)
        }
        else if( newValue < capacity ) {
            capacityDec( capacity-newValue )
        }
    }

    /**
     * @return the current capacity value by the number of slots specified
     */
    public int getCapacity() { capacity }

    /**
     * Increment the monitor capacity b
     *
     * @param slots The number of slots to add
     * @return The new capacity value
     */
    protected int capacityInc( int slots = 1 ) {
        def result = 0
        tasksQueueLock.withLock {
            result = capacity += slots
            notFull.signal()
        }
        log.debug "Monitor current capacity: $result (after inc: $slots)"
        return result
    }

    /**
     * Decrement the monitor capacity by the number of slots specified
     *
     * @param slots
     * @return The new capacity value
     */
    protected int capacityDec( int slots = 1 ) {

        def result = 0
        tasksQueueLock.withLock {
            result = capacity -= slots
        }

        log.debug "Monitor current capacity: $result (after dec: $slots)"
        return result
    }

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
        pollingQueue.size() < capacity
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
        pollingQueue.add(handler)
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
        pollingQueue.remove(handler)
    }

    /**
     * Schedule a new task for execution
     *
     * @param handler
     *      A {@link TaskHandler} representing the task to be submitted for execution
     */
    @Override
    void schedule(TaskHandler handler) {
        //
        // This guarantee that the 'pollingQueue' does not contain
        // more entries than the specified 'capacity'
        //
        boolean done = tasksQueueLock.withLock(true) {

            while ( !canSubmit(handler) )
                notFull.await();

            if( !session.isTerminated() && !session.isCancelled() ) {
                submit(handler)
                return true
            }

            return false
        }

        if( done )
            session.notifyTaskSubmit(handler)
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

        tasksQueueLock.withLock {
            if( remove(handler) ) {
                notFull.signal()
                return true
            }

            return false
        }
    }


    /**
     * Launch the monitoring thread
     *
     * @return
     *      The monitor object itself
     */
    @Override
    TaskMonitor start() {
        log.debug ">>> barrier register (monitor: ${this.name})"
        session.barrier.register(this)

        // creates the lock and condition
        this.tasksQueueLock = new ReentrantLock()
        this.notFull = tasksQueueLock.newCondition()

        this.taskCompleteLock = new ReentrantLock()
        this.taskComplete = taskCompleteLock.newCondition()

        // remove pending tasks on termination
        session.onShutdown { this.cleanup() }

        // launch the thread polling the queue
        Thread.start {
            try {
                pollLoop()
            }
            finally {
                log.debug "<<< barrier arrives (monitor: ${this.name})"
                session.barrier.arrive(this)
            }
        }

        return this
    }

    /**
     * Implements the polling strategy
     */
    protected void pollLoop() {

        while( true ) {
            long time = System.currentTimeMillis()
            log.trace "Scheduler queue size: ${pollingQueue.size()}"

            // check all running tasks for termination
            checkAllTasks()

            if( (session.isTerminated() && pollingQueue.size() == 0) || session.isAborted() ) {
                break
            }

            await(time)

            if( session.isAborted() ) {
                break
            }

            // dump this line every two minutes
            dumpInterval.throttle(true) {
                log.debug "!! executor $name > tasks to be completed: ${pollingQueue.size()} -- first: ${pollingQueue.peek()}"
            }
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


    /**
     * Check and update the status of queued tasks
     */
    protected void checkAllTasks() {

        for( TaskHandler handler : pollingQueue ) {
            try {
                checkTaskStatus(handler)
            }
            catch (Throwable error) {
                handleException(handler, error)
            }
        }

    }


    final protected void handleException( TaskHandler handler, Throwable error ) {
        def fault = null
        try {
            fault = handler.task.processor.resumeOrDie(handler?.task, error)
        }
        finally {
            // abort the session if a task task was returned
            if (fault instanceof TaskFault) {
                session.fault(fault)
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
        if( handler.checkIfCompleted()  ) {
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
                session.fault(fault)
            }
        }

    }

    /**
     * Kill all pending jobs when current execution session is aborted
     */
    protected void cleanup() {
        if( !pollingQueue.size() ) return
        log.warn "Killing pending tasks (${pollingQueue.size()})"

        def batch = new BatchCleanup()
        while( pollingQueue.size() ) {

            TaskHandler handler = pollingQueue.poll()
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
}

