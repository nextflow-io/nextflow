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

import groovy.util.logging.Slf4j
import nextflow.Session
/**
 * Monitor tasks execution for completion notifying their results
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class TaskDispatcher {

    final private Session session

    TaskDispatcher( Session session ) {

        this.session = session

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
        final handler = task.processor.executor.createTaskHandler(task)
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
     * Notify task start event
     */
    public void notifyStarted(TaskHandler handler) {

    }

    /**
     * Notify task termination event
     *
     * @param handler
     */
    public void notifyTerminated(TaskHandler handler) {

        // finalize the tasks execution
        handler.task.processor.finalizeTask(handler.task)
    }

    /**
     * Notify a task failure
     *
     * @param handler
     * @param e
     */
    public void notifyError(TaskHandler handler, Throwable e ) {
        handler.task.processor.handleException(e, handler.task)
    }



}

