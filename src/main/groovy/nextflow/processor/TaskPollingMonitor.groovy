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

    /**
     * Initialise the monitor, creating the queue which will hold the tasks
     *
     * @param session
     * @param execName
     * @param defPollInterval
     */
    TaskPollingMonitor( Session session, int capacity, Duration pollInterval, Duration dumpInterval = Duration.of('5min') ) {
        assert session, "Session object cannot be null"

        this.session = session
        this.dispatcher = session.dispatcher
        this.pollIntervalMillis = pollInterval.toMillis()
        this.dumpInterval = dumpInterval
        this.capacity = capacity

        this.pollingQueue = new ConcurrentLinkedQueue<>()

    }


    static TaskPollingMonitor create( Session session, String name, int defQueueSize, Duration defPollInterval ) {
        assert session
        assert name
        final capacity = session.getQueueSize(name, defQueueSize)
        final pollInterval = session.getPollInterval(name, defPollInterval)
        final dumpInterval = session.getMonitorDumpInterval(name)

        log.debug "Creating task monitor for executor '$name' > capacity: $capacity; pollInterval: $pollInterval; dumpInterval: $dumpInterval "
        new TaskPollingMonitor(session, capacity, pollInterval, dumpInterval)
    }

    static TaskPollingMonitor create( Session session, String name, Duration defPollInterval ) {
        assert session
        assert name

        final pollInterval = session.getPollInterval(name, defPollInterval)
        final dumpInterval = session.getMonitorDumpInterval(name)

        log.debug "Creating task monitor for executor '$name' > pollInterval: $pollInterval; dumpInterval: $dumpInterval "
        new TaskPollingMonitor(session, 0, pollInterval, dumpInterval)
    }

    protected Queue<TaskHandler> getPollingQueue() { pollingQueue }

    /**
     * Set the number of tasks the monitor will handle in parallel.
     * This value must be changed by using a mutex lock
     *
     * @param value
     */
    protected void setCapacity( int value ) {
        this.capacity = value
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

        mutex.withLock {
            def result = capacity += slots
            notFull.signal()
            return result
        }

    }

    /**
     * Decrement the monitor capacity by the number of slots specified
     *
     * @param slots
     * @return The new capacity value
     */
    protected int capacityDec( int slots = 1 ) {

        mutex.withLock {
            return capacity -= slots
        }
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
        mutex.withLock(true) {

            while ( pollingQueue.size() >= capacity )
                notFull.await();

            handler.submit()
            // the entry is appended to the 'submitQueue', not directly to 'pollingQueue'

            pollingQueue.add(handler)
        }

        dispatcher.notifySubmitted(handler)
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
    TaskHandler getTaskHandlerBy( final taskId ) {
        pollingQueue.find { handler -> handler.task.id == taskId }
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

            // check all running tasks for termination
            checkAll(pollingQueue)

            if( session.isTerminated() && pollingQueue.size() == 0 ) {
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
    void signalComplete() {
        mutex.withLock {
            taskComplete.signal()
        }
    }


    /**
     * @param collection The collections of tasks to check
     */
    protected void checkAll( Collection<TaskHandler> collection ) {

        collection.each { handler ->

            try {
                checkTaskStatus(handler)
            }
            catch( Exception e ) {
                dispatcher.notifyError(e, handler)
            }

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
            dispatcher.notifyStarted(handler)
        }

        // check if it is terminated
        if( handler.checkIfCompleted()  ) {
            // since completed *remove* the task from the processing queue
            drop(handler)
            // finalize the tasks execution
            dispatcher.notifyTerminated(handler)
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

