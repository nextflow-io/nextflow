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

    final BlockingQueue<TaskRun> queue

    final private Session session

    TaskScheduler( Session session ) {

        this.session = session
        this.queue = new ArrayBlockingQueue<TaskRun>( getQueueSize(session.config) )

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
     * Add the specify task instance to the queue of launched tasks
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
        queue.put(task)

        try {
            task.processor.launchTask(task)
        }
        catch( Exception error ) {
            queue.remove(task)
            throw error
        }
    }


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


    private void checkAll( Collection<TaskRun> collection ) {
        collection.each { TaskRun task ->
            checkTask(task)
        }
    }

    private void checkTask( TaskRun task ) {
        assert task

        task.processor.with {

            // check if it is started
            def started = checkTaskStarted(task)

            // check if it is terminated
            if( started && checkTaskCompletion(task)  ) {

                // since completed *remove* the task from the processing queue
                log.debug "Removing task from running queue: ${task}"
                queue.remove(task)

                // finalize the tasks execution
                finalizeTask(task)
            }
        }
    }


    private void killOnExit()  {

        Runtime.addShutdownHook {
            while( queue.size() ) {
                TaskRun task = queue.poll()
                log.warn "Killing pending task: $task"
                task.handler?.kill()
            }
        }

    }

}

