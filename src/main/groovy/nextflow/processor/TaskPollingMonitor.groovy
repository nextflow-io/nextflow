/*
 * Copyright (c) 2013-2015, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2015, Paolo Di Tommaso and the respective authors.
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

import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.util.Duration
/**
 *
 * Monitors the queued tasks waiting for their termination
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@Slf4j
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
    private Lock mutex

    /**
     * A condition that signal when at least a complete task is available
     */
    private Condition taskComplete

    private Condition notFull

    private Queue<TaskHandler> pollingQueue

    private int capacity

    private Queue<Closure> eventsQueue

    private Queue<Closure> listenersQueue

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
    TaskPollingMonitor( Map params ) {
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
        this.eventsQueue = new ConcurrentLinkedQueue<>()
        this.listenersQueue = new ConcurrentLinkedQueue<>()
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
        mutex.withLock {
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
        mutex.withLock {
            result = capacity -= slots
        }

        log.debug "Monitor current capacity: $result (after dec: $slots)"
        return result
    }


    /**
     * Add a new task handler to queue tasks queue
     *
     * @param handler
     */
    @Override
    void put(TaskHandler handler) {

        //
        // This guarantee that the 'pollingQueue' does not contain
        // more entries than the specified 'capacity'
        //
        boolean done = false
        mutex.withLock(true) {

            while ( pollingQueue.size() >= capacity )
                notFull.await();

            if( !session.isTerminated()) {
                // submit the job execution -- throws a ProcessException when submit operation fail
                handler.submit()
                // note: add the 'handler' into the polling queue *after* the submit operation,
                // this guarantees that in the queue are only jobs successfully submitted
                pollingQueue.add(handler)
                done = true
            }

        }

        if( done )
            dispatcher.notifySubmit(handler)
    }

    /**
     * Remove a task from the processing tasks queue
     *
     * @param handler A {@code HzTaskHandler} instance
     * @return
     */
    @Override
    boolean drop(TaskHandler handler) {
        if( !handler ) {
            log.debug "Unknown task handler to drop: $handler"
            return false
        }

        log.trace "Dropping task handler: $handler"
        mutex.withLock {
            if( pollingQueue.remove(handler) ) {
                notFull.signal()
                return true
            }

            return false
        }
    }

    /**
     * Get a queued {@code TaskHandler} by the {@code TaskRun#id}
     *
     * @param taskId
     * @return
     */
    TaskHandler getTaskHandlerById( final taskId ) {
        pollingQueue.find { handler -> handler.task.id == taskId }
    }

    /**
     * Find all {@code TaskHandler} in the specified status
     *
     * @param status The status of the required task handlers
     * @return The list of {@code TaskHandler} win the specified {@code TaskHandler.Status}, or an empty list if not tashs are found.
     */
    TaskHandler getTaskHandlerByStatus( final TaskStatus status ) {
        pollingQueue.find { handler -> handler.status == status }
    }

    /**
     * Find all the task handlers that matches the specified closure
     * @param closure
     * @return
     */
    List<TaskHandler> findTaskHandlers( Closure closure = null ) {
        closure ? pollingQueue.findAll(closure) : pollingQueue.findAll()
    }


    /**
     * Launch the monitoring thread
     *
     * @return The monitor object itself
     */
    def TaskMonitor start() {
        log.debug ">>> phaser register (scheduler)"
        session.phaser.register()

        // creates the lock and condition
        this.mutex = new ReentrantLock()
        this.notFull = mutex.newCondition()
        this.taskComplete = mutex.newCondition()

        // remove pending tasks on termination
        session.onShutdown { this.cleanup() }

        // launch the thread polling the queue
        Thread.start {
            try {
                pollLoop()
            }
            finally {
                log.debug "<<< phaser de-register (scheduler)"
                session.phaser.arriveAndDeregister()
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

            // process all scheduled events
            processEvents()

            // check all running tasks for termination
            checkAllTasks()

            processListeners()

            if( (session.isTerminated() && pollingQueue.size() == 0) || session.isAborted() ) {
                break
            }

            await(time)

            // dump this line every two minutes
            dumpInterval.throttle(true) {
                log.debug "!! executor $name > tasks to be completed: ${pollingQueue.size()} -- first: ${pollingQueue.peek()}"
            }
        }

    }

    /**
     * Await for one or more tasks to be processed
     *
     * @param time The wait timeout in millis
     */
    protected void await( long time ) {
        def delta = this.pollIntervalMillis - (System.currentTimeMillis() - time)
        if( delta <= 0 )
            return

        mutex.withLock {
            taskComplete.await( delta, TimeUnit.MILLISECONDS )
        }
    }

    /**
     * Signal that a task has been completed
     */
    void signal() {
        mutex.withLock {
            taskComplete.signal()
        }
    }


    /**
     * Check and update the queue tasks status
     *
     * @param collection The collections of tasks to check
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

    /**
     * Applies and consumes all the queued events
     */
    final protected processEvents() {

        Iterator<Closure> itr = eventsQueue.iterator()
        while( itr.hasNext() ) {
            Closure event = itr.next()
            try {
                event.call()
            }
            finally {
                itr.remove()
            }
        }

    }

    /**
     * Process all registered monitor listeners
     *
     * See #register
     */
    final protected processListeners() {

        for( Closure closure : listenersQueue ) {
            try {
                closure.call(this)
            }
            catch( Exception e ) {
                log.error "Error processing monitor listener", e
            }
        }

    }

    /**
     * Add a new listener to the task monitor
     *
     * @param listener
     * @return
     */
    final void register( Closure listener ) {
        if( !listener ) return
        listenersQueue << listener
    }


    final protected void handleException( TaskHandler handler, Throwable error ) {
        try {
            handler.task.processor.resumeOrDie(handler?.task, error)
        }
        finally {
            dispatcher.notifyError(handler, error)
        }
    }


    /**
     * Check the status of the given task
     *
     * @param handler The {@code TaskHandler} instance of the task to check
     */
    protected void checkTaskStatus( TaskHandler handler ) {
        assert handler

        // check if it is started
        if( handler.checkIfRunning() ) {
            dispatcher.notifyStart(handler)
        }

        // check if it is terminated
        if( handler.checkIfCompleted()  ) {
            // since completed *remove* the task from the processing queue
            drop(handler)

            // finalize the tasks execution
            handler.task.processor.finalizeTask(handler.task)
            // trigger the count down latch when it is a blocking task
            handler.latch?.countDown()

            // finalize the tasks execution
            dispatcher.notifyComplete(handler)
        }

    }


    /**
     * Kill all pending jobs when aborting
     */

    protected void cleanup () {
        if( !pollingQueue.size() ) return
        log.warn "Killing pending processes (${pollingQueue.size()})"

        int c = 0
        while( pollingQueue.size() ) {
            TaskHandler handler = pollingQueue.poll()
            try {
                log.debug "Killing process (${++c}): $handler"
                handler.kill()
            }
            catch( Throwable e ) {
                log.debug "Failed killing pending process: ${handler} -- cause: ${e.message}"
            }
        }
    }
}

