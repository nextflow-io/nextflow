package nextflow.processor

import java.util.concurrent.ArrayBlockingQueue
import java.util.concurrent.BlockingQueue

import groovy.util.logging.Slf4j
import nextflow.Session
/**
 *
 * Maintains a queue of tasks for which manage status transition
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@Slf4j
class TaskPollingQueue implements TaskQueueHolder {

    final BlockingQueue<TaskHandler> queue

    final Session session

    final TaskDispatcher dispatcher

    final long pollInterval


    TaskPollingQueue( Session session, long pollInterval ) {
        assert session
        assert pollInterval

        this.session = session
        this.dispatcher = session.dispatcher
        this.pollInterval = pollInterval
        this.queue = new ArrayBlockingQueue<TaskHandler>( getQueueSize(session.config) )

        killOnExit()
    }



    private int getQueueSize( Map config ) {
        // creating the running tasks queue
        def size
        if( config.queueSize ) {
            size = config.queueSize
            log.debug "Processor runnable queue size: $size"
        }
        else {
            size = Math.max( Runtime.getRuntime().availableProcessors()-1, 1 )
            log.warn "Undefined processor runnable queue size -- fallback on num of available processors-1: $size"
        }
        return size
    }


    def void start() {
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

    }

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


    @Override
    def void put( TaskHandler handler ) {
        queue.put( handler )
    }

    @Override
    def boolean remove( TaskHandler handler ) {
        queue.remove(handler)
    }


    private void checkAll( Collection<TaskHandler> collection ) {
        collection.each { TaskHandler handler ->
            try {
                checkTask(handler)
            }
            catch( Exception e ) {
                dispatcher.notifyError(handler, e)
            }
        }
    }

    /**
     * Check the status of the given task
     *
     * @param handler The {@code TaskHandler} instance of the task to check
     */
    private void checkTask( TaskHandler handler ) {
        assert handler

        // check if it is started
        def started = handler.isStarted()
        if( !started && handler.checkIfStarted() ) {
            dispatcher.notifyStarted(handler)
            started = true
        }

        // check if it is terminated
        if( started && handler.checkIfTerminated()  ) {
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
