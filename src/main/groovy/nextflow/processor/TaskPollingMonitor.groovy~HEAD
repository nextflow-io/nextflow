package nextflow.processor

import java.util.concurrent.ArrayBlockingQueue
import java.util.concurrent.BlockingQueue

import groovy.util.logging.Slf4j
import nextflow.Session
/**
 *
 * Monitors the queued tasks waiting for their termination
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@Slf4j
class TaskPollingMonitor implements TaskMonitor {


    /**
     * The submitted tasks queue. It is implemented using a blocking queue,
     * so that no more than a fixed amount of tasks can be submitted concurrently
     */
    final BlockingQueue<TaskHandler> queue

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
    final long pollInterval

    /**
     * Initialise the monitor, creatign the queue which will hold the tasks
     *
     * @param session
     * @param queueSize
     * @param pollInterval
     */
    TaskPollingMonitor( Session session, int queueSize, long pollInterval ) {
        assert session , "Session object cannot be null"
        assert queueSize, "Executor queue size cannot be zero"
        assert pollInterval, "Executor polling interval cannot be zero"

        this.session = session
        this.dispatcher = session.dispatcher
        this.pollInterval = pollInterval
        this.queue = new ArrayBlockingQueue<TaskHandler>(queueSize)

        killOnExit()
    }

    /**
     * Launch the monitoring thread
     *
     * @return The monitor object itself
     */
    def TaskPollingMonitor start() {
        log.debug ">>> phaser register (scheduler)"
        session.phaser.register()

        Thread.start {
            try {
                pollLoop()
            }
            finally {
                log.debug "<<< phaser deregister (scheduler)"
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
            log.trace "Scheduler queue size: ${queue.size()}"

            // check all running tasks for termination
            checkAll(queue)

            if( session.isTerminated() && queue.size() == 0 ) {
                break
            }

            def delta = pollInterval - (System.currentTimeMillis() - time)
            if( delta>0 ) { sleep(delta) }

        }
    }

    /**
     * Add an entry to the queue of monitored tasks
     *
     * @param handler The {@code TaskHandler} instance for the task
     */
    @Override
    def void put( TaskHandler handler ) {
        queue.put( handler )
    }

    /**
     * Remove an entry from the queue of monitored tasks
     *
     * @param handler The {@code TaskHandler} instance for the task
     */
    @Override
    def boolean remove( TaskHandler handler ) {
        queue.remove(handler)
    }

    /**
     * @param collection The collections of tasks to check
     */
    private void checkAll( Collection<TaskHandler> collection ) {
        collection.each { TaskHandler handler ->
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
    private void checkTaskStatus( TaskHandler handler ) {
        assert handler

        // check if it is started
        if( handler.checkIfRunning() ) {
            dispatcher.notifyStarted(handler)
        }

        // check if it is terminated
        if( handler.checkIfCompleted()  ) {
            // since completed *remove* the task from the processing queue
            queue.remove(handler)
            // finalize the tasks execution
            dispatcher.notifyTerminated(handler)
        }

    }


    /**
     * Kill all pending jobs when aborting
     */
    private void killOnExit()  {

        Runtime.addShutdownHook {
            while( queue.size() ) {
                TaskHandler handler = queue.poll()
                log.warn "Killing pending task: ${handler.task.name}"
                handler.kill()
            }
        }

    }
}
