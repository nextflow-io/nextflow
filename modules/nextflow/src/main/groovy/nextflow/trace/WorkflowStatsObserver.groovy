/*
 * Copyright 2013-2026, Seqera Labs
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

package nextflow.trace

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import groovyx.gpars.agent.Agent
import nextflow.Session
import nextflow.processor.TaskProcessor
import nextflow.trace.event.TaskEvent
import nextflow.util.Threads

/**
 * Collect workflow execution progress and statistics
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class WorkflowStatsObserver implements TraceObserverV2 {

    // updated by the agent, and copied under its monitor
    private final WorkflowStats data = new WorkflowStats()

    private Session session

    private Agent<WorkflowStats> agent

    WorkflowStatsObserver(Session session) {
        this.session = session
        this.agent = new Agent<>(data)
    }

    /**
     * Runs the action on the agent. An error aborts the run.
     */
    private void send(Closure action) {
        agent.send {
            try {
                synchronized (data) {
                    action.call()
                }
            }
            catch( Throwable e ) {
                log.debug "Unexpected error while updating the workflow stats | ${e.message ?: e}"
                // abort on another thread, since the completion handlers wait for the agent
                Threads.start { session.abort(e) }
            }
        }
    }

    @Override
    void onFlowComplete() {
        log.debug "Workflow completed > ${getQuickStats()}"
    }

    @Override
    void onProcessCreate(TaskProcessor process) {
        log.trace "== event create pid=${process.id}"
        send { data.markCreated(process) }
    }

    @Override
    void onProcessTerminate(TaskProcessor process) {
        log.trace "== event terminated pid=${process.id}"
        send { data.markTerminated(process) }
    }

    @Override
    void onTaskPending(TaskEvent event) {
        log.trace "== event pending pid=${event.handler.task.processor.id}; status=$event.handler.status"
        send { data.markPending(event.handler.task.processor) }
    }

    @Override
    void onTaskSubmit(TaskEvent event) {
        log.trace "== event submit pid=${event.handler.task.processor.id}; status=$event.handler.status"
        send { data.markSubmitted(event.handler.task) }
    }

    @Override
    void onTaskStart(TaskEvent event) {
        log.trace "== event start pid=${event.handler.task.processor.id}; status=$event.handler.status"
        send { data.markRunning(event.handler.task) }
    }

    @Override
    void onTaskComplete(TaskEvent event) {
        log.trace "== event complete pid=${event.handler.task.processor.id}; status=$event.handler.status"
        send { data.markCompleted(event.handler.task, event.trace, event.handler.status) }
    }

    @Override
    void onTaskCached(TaskEvent event) {
        log.trace "== event cached pid=${event.handler.task.processor.id}; status=$event.handler.status"
        send { data.markCached(event.handler.task, event.trace) }
    }

    /**
     * A copy of the stats, once the events notified before are applied,
     * or without waiting for them when the caller is interrupted or
     * runs on the agent thread
     */
    WorkflowStats getStats() {
        try {
            agent.await()
        }
        catch( InterruptedException e ) {
            log.warn "Got an interrupted exception while waiting for the workflow stats | ${e}"
            Thread.currentThread().interrupt()
        }
        catch( IllegalStateException e ) {
            log.debug "Unable to await the workflow stats from the agent thread | ${e.message}"
        }
        return getQuickStats()
    }

    /**
     * A copy of the stats, without waiting for the pending events
     */
    WorkflowStats getQuickStats() {
        synchronized (data) {
            return data.clone()
        }
    }

    boolean hasProgressRecords() {
        return data.getProgressLength()
    }

    long getChangeTimestamp() {
        return data.changeTimestamp
    }
}
