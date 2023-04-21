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

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.executor.ArrayTaskHandler
import nextflow.util.Duration

/**
 * Extension of the polling monitor for array jobs.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class ArrayTaskPollingMonitor extends TaskPollingMonitor {

    protected ArrayTaskPollingMonitor(Map params) {
        super(params)
    }

    static ArrayTaskPollingMonitor create( Session session, String name, int defQueueSize, Duration defPollInterval ) {
        final capacity = session.getQueueSize(name, defQueueSize)
        final pollInterval = session.getPollInterval(name, defPollInterval)
        final dumpInterval = session.getMonitorDumpInterval(name)

        log.debug "Creating task monitor for executor '$name' > capacity: $capacity; pollInterval: $pollInterval; dumpInterval: $dumpInterval "
        new ArrayTaskPollingMonitor(name: name, session: session, capacity: capacity, pollInterval: pollInterval, dumpInterval: dumpInterval)
    }

    protected ArrayTaskHandler toArrayHandler(TaskHandler handler) {
        if( handler !instanceof ArrayTaskHandler )
            throw new IllegalStateException()

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
            handler.array.each { h ->
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
            handler.array.each { h ->
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

