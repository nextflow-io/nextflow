/*
 * Copyright (c) 2012, the authors.
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
import java.util.concurrent.ArrayBlockingQueue
import java.util.concurrent.BlockingQueue

import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.util.Duration
/**
 * Monitor tasks execution for completion notifying their results
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class TaskScheduler {

    final BlockingQueue<TaskHandler> queue

    final private Session session

    TaskScheduler( Session session ) {

        this.session = session
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
            size = Runtime.getRuntime().availableProcessors()
            log.warn "Undefined processor runnable queue size -- fallback on num of available processors: $size"
        }
        return size
    }

    /**
     * Submit the specified task for execution to the underlying system
     * and add it to the queue of tasks to be monitored.
     *
     * @param task A {@code TaskRun} instance
     */
    void submit( TaskRun task ) {
        log.debug "Scheduling task: ${task}"
        if( session.isTerminated() ) {
            new IllegalStateException("Session terminated - Cannot add task to execution queue: ${task}")
        }

        /*
         * Add the task to the queue for processing
         * Note: queue is implemented as a fixed size blocking queue, when
         * there's not space *put* operation will block until, some other tasks finish
         */
        def handler = task.processor.executor.createTaskHandler(task)
        queue.put(handler)

        try {
            handler.submit()
        }
        catch( Exception error ) {
            queue.remove(handler)
            throw error
        }
    }


    /**
     * Launch the scheduler which will monitor all
     */
    public TaskScheduler start () {
        log.debug ">>> phaser register (scheduler)"
        session.phaser.register()

        Thread.start('Tasks scheduler') {

            log.debug "Starting tasks scheduler"

            def MAX_WAIT = Duration.create('1s').toMillis()
            while( true ) {
                long time = System.currentTimeMillis()
                log.trace "Scheduler queue size: ${queue.size()}"

                // check all running tasks for termination
                checkAll(queue)

                if( session.isTerminated() && queue.size() == 0 ) {
                    break
                }

                def delta = MAX_WAIT - (System.currentTimeMillis() - time)
                if( delta>0 ) { sleep(delta) }

            }

            log.debug "<<< phaser deregister (scheduler)"
            session.phaser.arriveAndDeregister()
            log.debug "Terminating tasks scheduler"

        }

        return this
    }


    private void checkAll( Collection<TaskHandler> collection ) {
        collection.each { TaskHandler handler ->
            checkTask(handler)
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
        def started = handler.checkIfRunning()

        // check if it is terminated
        if( started && handler.checkIfTerminated()  ) {

            // since completed *remove* the task from the processing queue
            log.debug "Removing task from running queue: ${handler}"
            queue.remove(handler)

            // finalize the tasks execution
            handler.task.processor.finalizeTask(handler.task)
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

