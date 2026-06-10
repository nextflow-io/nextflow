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
import nextflow.Session
import nextflow.processor.TaskProcessor
import nextflow.trace.event.TaskEvent
import nextflow.util.SimpleAgent

/**
 * Collect workflow execution progress and statistics
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class WorkflowStatsObserver implements TraceObserverV2 {

    private WorkflowStats data = new WorkflowStats()

    private Session session

    private SimpleAgent<WorkflowStats> agent

    WorkflowStatsObserver(Session session) {
        this.session = session
        this.agent = new SimpleAgent(data).onError { err -> session.abort(err) }
    }

    @Override
    void onFlowComplete() {
        log.debug "Workflow completed > $data"
    }

    @Override
    void onProcessCreate(TaskProcessor process) {
        log.trace "== event create pid=${process.id}"
        agent.send { data.markCreated(process) }
    }

    @Override
    void onProcessTerminate(TaskProcessor process) {
        log.trace "== event terminated pid=${process.id}"
        agent.send { data.markTerminated(process) }
    }

    @Override
    void onTaskPending(TaskEvent event) {
        log.trace "== event pending pid=${event.handler.task.processor.id}; status=$event.handler.status"
        agent.send { data.markPending(event.handler.task.processor) }
    }

    @Override
    void onTaskSubmit(TaskEvent event) {
        log.trace "== event submit pid=${event.handler.task.processor.id}; status=$event.handler.status"
        agent.send { data.markSubmitted(event.handler.task) }
    }

    @Override
    void onTaskStart(TaskEvent event) {
        log.trace "== event start pid=${event.handler.task.processor.id}; status=$event.handler.status"
        agent.send { data.markRunning(event.handler.task) }
    }

    @Override
    void onTaskComplete(TaskEvent event) {
        log.trace "== event complete pid=${event.handler.task.processor.id}; status=$event.handler.status"
        agent.send { data.markCompleted(event.handler.task, event.trace, event.handler.status) }
    }

    @Override
    void onTaskCached(TaskEvent event) {
        log.trace "== event cached pid=${event.handler.task.processor.id}; status=$event.handler.status"
        agent.send { data.markCached(event.handler.task, event.trace) }
    }

    WorkflowStats getStats() {
        return agent.getValue()
    }

    WorkflowStats getQuickStats() {
        return agent.getQuickValue()
    }

    boolean hasProgressRecords() {
        return data.getProgressLength()
    }

    long getChangeTimestamp() {
        return data.changeTimestamp
    }
}
