/*
 * Copyright 2013-2023, Seqera Labs
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

import java.util.concurrent.ConcurrentHashMap
import java.util.concurrent.ConcurrentLinkedQueue

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.executor.ArrayTaskAware
import nextflow.executor.ArrayTaskHandler
import nextflow.executor.Executor
import nextflow.trace.TraceObserver
import nextflow.util.Duration

/**
 * Polling monitor that submits tasks in batches.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class ArrayTaskPollingMonitor extends TaskPollingMonitor implements TraceObserver {

    private Executor executor

    private int arraySize

    private Map<String,Queue<TaskHandler>> arrayQueues = new ConcurrentHashMap<>()

    private Map<String,Boolean> closedProcesses = new ConcurrentHashMap<>()

    protected ArrayTaskPollingMonitor(Map params) {
        super(params)
        this.executor = params.executor as Executor
        this.arraySize = params.arraySize as int

        session.registerObserver(this)
    }

    static ArrayTaskPollingMonitor create( Session session, Executor executor, int defQueueSize, Duration defPollInterval, int defArraySize ) {
        final name = 'array'
        final capacity = session.getQueueSize(name, defQueueSize)
        final pollInterval = session.getPollInterval(name, defPollInterval)
        final dumpInterval = session.getMonitorDumpInterval(name)
        final arraySize = session.getExecConfigProp(name, 'arraySize', defArraySize) as Integer

        log.debug "Creating array task monitor for exector '$executor.name' > capacity: $capacity; pollInterval: $pollInterval; dumpInterval: $dumpInterval; arraySize: $arraySize"

        new ArrayTaskPollingMonitor(
            name: name,
            session: session,
            executor: executor,
            capacity: capacity,
            pollInterval: pollInterval,
            dumpInterval: dumpInterval,
            arraySize: arraySize )
    }

    /**
     * Add scheduled tasks to a queue, and schedule an array job when
     * the queue reaches the desired size.
     *
     * @param handler
     */
    @Override
    synchronized void schedule(TaskHandler handler) {
        final process = handler.task.processor.name

        // schedule task directly if process has already closed
        if( process in closedProcesses ) {
            executor.monitor.schedule(handler)
            return
        }

        // initialize array queue
        if( process !in arrayQueues )
            arrayQueues[process] = new ConcurrentLinkedQueue<>()

        // add task to the array queue
        final queue = arrayQueues[process]
        queue.add(handler)

        // schedule array job when a batch is ready
        if( queue.size() >= arraySize ) {
            log.debug "Submitting array job for process '${process}'"
            schedule0(queue, arraySize)
        }
    }

    /**
     * Submit any remaining tasks as an array job when a process is closed.
     *
     * @param process
     */
    @Override
    void onProcessClose(String process) {
        final queue = arrayQueues[process]

        if( queue != null && queue.size() > 0 ) {
            log.debug "Submitting remainder array job for process '${process}'"
            schedule0(queue, queue.size())
        }

        closedProcesses[process] = true
    }

    synchronized protected void schedule0( Queue<TaskHandler> queue, int size ) {
        // remove tasks from the queue
        def array = new ArrayList<TaskHandler>()
        def iter = queue.iterator()

        while( iter.hasNext() && array.size() < size ) {
            array << iter.next()
            iter.remove()
        }

        if( array.size() < size )
            log.warn "Array job expected ${size} tasks but received only ${array.size()}"

        // create array task handler and schedule it
        super.schedule(((ArrayTaskAware)executor).createArrayTaskHandler(array))
    }

    protected ArrayTaskHandler toArrayHandler(TaskHandler handler) {
        assert handler instanceof ArrayTaskHandler
        return (ArrayTaskHandler)handler
    }

    @Override
    protected void handleException(TaskHandler handler, Throwable error) {
        handler = toArrayHandler(handler)

        def fault = null
        def faultHandler = null
        try {
            // remove the array task from the processing queue
            if( evict(handler) )
                handler.decProcessForks()

            // attempt to retry each task in the array
            for( TaskHandler h : handler.array ) {
                def fault0 = h.task.processor.resumeOrDie(h.task, error)
                if( fault == null && fault0 ) {
                    fault = fault0
                    faultHandler = h
                }
            }

            log.trace "Task fault (1): $fault"
        }
        finally {
            // abort the session if a task fault was returned
            if( fault instanceof TaskFault )
                session.fault(fault, handler)
        }
    }

    @Override
    protected void checkTaskStatus( TaskHandler handler ) {
        handler = toArrayHandler(handler)

        // check if the array task is started
        if( handler.checkIfRunning() ) {
            log.trace "Task started > $handler"
            notifyTaskStart(handler)
        }

        // check if the array task is completed
        if( handler.checkIfCompleted() ) {
            log.debug "Task completed > $handler"

            // decrement forks count
            handler.decProcessForks()

            // remove the array task from the processing queue
            evict(handler)

            // attempt to finalize each task in the array
            def fault = null
            def faultHandler = null
            for( TaskHandler h : handler.array ) {
                final fault0 = h.task.processor.finalizeTask(h.task)
                if( fault == null && fault0 ) {
                    fault = fault0
                    faultHandler = h
                }
            }

            // notify task completion
            notifyTaskComplete(handler)

            // abort the session if a task fault was returned
            if( fault instanceof TaskFault )
                session.fault(fault, faultHandler)
        }
    }

    @Override
    protected void notifyTaskSubmit(TaskHandler handler) {
        toArrayHandler(handler).array.each { h -> session.notifyTaskSubmit(h) }
    }

    @Override
    protected void notifyTaskPending(TaskHandler handler) {
        toArrayHandler(handler).array.each { h -> session.notifyTaskPending(h) }
    }

    @Override
    protected void notifyTaskStart(TaskHandler handler) {
        toArrayHandler(handler).array.each { h -> session.notifyTaskStart(h) }
    }

    @Override
    protected void notifyTaskComplete(TaskHandler handler) {
        toArrayHandler(handler).array.each { h -> session.notifyTaskComplete(h) }
    }

}
