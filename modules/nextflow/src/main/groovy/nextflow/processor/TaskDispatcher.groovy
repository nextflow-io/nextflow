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
import java.util.concurrent.CountDownLatch

import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.executor.Executor
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
    void submit( TaskRun task, boolean awaitTermination ) {
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
        monitor.schedule(handler)
        if( handler && handler.latch ) {
            log.trace "Process ${task} > blocking"
            handler.latch.await()
            log.trace "Process ${task} > complete"
        }
    }

}

