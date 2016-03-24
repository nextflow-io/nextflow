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
import java.util.concurrent.CountDownLatch

import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.executor.Executor
import nextflow.trace.TraceObserver
/**
 * Monitor tasks execution for completion notifying their results
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class TaskDispatcher {

    /**
     * The current session object
     */
    final private Session session

    /**
     * Map each executor class with its tasks monitor, in other words there's one {@code TaskMonitor}
     * instance for each type of executor
     */
    final private Map<Class<? extends Executor>, TaskMonitor> monitors = [:]

    private volatile boolean started

    /*
     * A list of {@link TraceObserver} instances to which notify process events
     */
    private List<TraceObserver> observers = []

    /**
     * Dispatcher constructor
     *
     * @param session
     */
    TaskDispatcher( Session session ) {
        this.session = session
    }


    /* ONLY FOR TEST PURPOSE */
    protected TaskDispatcher( ) {}


    /**
     *
     * Get an instance of {@code TaskMonitor} in the {@code #monitors} map,
     * when the map contains no value for the specified executor invoke the
     * closure in order to create a new {@code TaskMonitor} object,
     * add it to the map and return it.
     *
     * NOTE: this method is not thread-safe by design since it is invoked
     * only by the script interpreter.
     *
     * @param type
     * @param create
     * @return
     */
    TaskMonitor getOrCreateMonitor( Class<? extends Executor> type, Closure<TaskMonitor> create ) {

        if( monitors.containsKey(type) ) {
            return monitors.get(type)
        }

        def result = create.call()
        monitors.put(type,result)

        if( started ) {
            log.debug "Starting monitor: ${result.class.simpleName}"
            result.start()
        }

        return result
    }

    /**
     * Get the monitor associated to the executor type or instance
     *
     * @param executor
     * @return
     */
    TaskMonitor getMonitor( executor ) {
        assert executor

        def result = null
        if( executor instanceof Executor )
            result = monitors.get(executor.class)

        else if( executor instanceof Class )
            result = monitors.get(executor)

        if( !result )
            throw new IllegalStateException("Missing monitor for executor: $executor")

        return result
    }

    /**
     * Flag the dispatcher as started
     */
    void start() {
        log.debug "Dispatcher > start"
        started = true
    }

    void signal() {
        monitors.values().each { TaskMonitor monitor -> monitor.signal() }
    }


    /**
     * Submit the specified task for execution to the underlying system
     * and add it to the queue of tasks to be monitored.
     *
     * @param task A {@code TaskRun} instance
     */
    void submit( TaskRun task, boolean awaitTermination, String submitMessage ) {
        log.trace "Scheduling process: ${task}"

        if( session.isTerminated() ) {
            new IllegalStateException("Session terminated - Cannot add process to execution queue: ${task}")
        }

        final executor = task.processor.executor
        final monitor = getMonitor(executor)
        final handler = executor.createTaskHandler(task)

        // set a count down latch if the execution is blocking
        if( awaitTermination )
            handler.latch = new CountDownLatch(1)

        /*
         * Add the task to the queue for processing
         * Note: queue is implemented as a fixed size blocking queue, when
         * there's not space *put* operation will block until, some other tasks finish
         */
        monitor.put(handler)
        log.info "[${task.hashLog}] $submitMessage > ${task.name}"
        if( handler && handler.latch ) {
            log.trace "Process ${task} > blocking"
            handler.latch.await()
            log.trace "Process ${task} > complete"
        }
    }

    /**
     * Notifies that a task has been submitted
     */
    public void notifySubmit( TaskHandler handler ) {
        for( TraceObserver it : observers ) {
            try {
                it.onProcessSubmit(handler)
            }
            catch( Exception e ) {
                log.error(e.getMessage(), e)
            }
        }
    }

    /**
     * Notifies task start event
     */
    public void notifyStart( TaskHandler handler ) {
        for( TraceObserver it : observers ) {
            try {
                it.onProcessStart(handler)
            }
            catch( Exception e ) {
                log.error(e.getMessage(), e)
            }
        }
    }

    /**
     * Notifies task termination event
     *
     * @param handler
     */
    public void notifyComplete( TaskHandler handler ) {
        // notify the event to the observers
        for( TraceObserver it : observers ) {
            try {
                it.onProcessComplete(handler)
            }
            catch( Exception e ) {
                log.error(e.getMessage(), e)
            }
        }
    }

    /**
     * Notify a task failure
     *
     * @param handler
     * @param e
     */
    public void notifyError( TaskHandler handler, Throwable error ) {
        for( TraceObserver it : observers ) {
            try {
                it.onProcessError(handler,error)
            }
            catch( Exception e ) {
                log.error(e.getMessage(), e)
            }
        }
    }

    /**
     * Register a {@link TraceObserver} instance to which process event will be notified
     *
     * @param obj An object implement {@link TraceObserver} trait
     */
    def void register( TraceObserver obj ) {
        // just return if it's a null object
        // or it is already contained in the list
        if( !obj || observers.contains(obj))
            return

        // append this obj to the observer list
        observers << obj
    }

    /**
     * Un-register an observer object from the list of registered observers
     *
     * @param obj An object implement {@link TraceObserver} trait
     */
    def void unregister( TraceObserver obj ) {
        if( !obj ) return
        observers.remove(obj)
    }


}

