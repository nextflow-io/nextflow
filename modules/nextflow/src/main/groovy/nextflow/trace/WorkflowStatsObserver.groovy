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

package nextflow.trace

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.processor.TaskHandler
import nextflow.processor.TaskProcessor
import nextflow.util.SimpleAgent

/**
 * Collect workflow execution progress and statistics
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class WorkflowStatsObserver implements TraceObserver {

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
    void onProcessCreate(TaskProcessor process){
        log.trace "== event create pid=${process.id}"
        agent.send { data.markCreated(process) }
    }

    @Override
    void onProcessTerminate( TaskProcessor processor ) {
        log.trace "== event terminated pid=${processor.id}"
        agent.send { data.markTerminated(processor) }
    }

    @Override
    void onProcessPending(TaskHandler handler, TraceRecord trace){
        log.trace "== event pending pid=${handler.getTask().processor.id}; status=$handler.status"
        agent.send { data.markPending(handler.getTask().processor) }
    }

    @Override
    void onProcessSubmit(TaskHandler handler, TraceRecord trace){
        log.trace "== event submit pid=${handler.getTask().processor.id}; status=$handler.status"
        agent.send { data.markSubmitted(handler.getTask()) }
    }

    @Override
    void onProcessStart(TaskHandler handler, TraceRecord trace){
        log.trace "== event start pid=${handler.getTask().processor.id}; status=$handler.status"
        agent.send { data.markRunning(handler.getTask()) }
    }

    @Override
    void onProcessComplete(TaskHandler handler, TraceRecord trace) {
        log.trace "== event complete pid=${handler.getTask().processor.id}; status=$handler.status"
        agent.send { data.markCompleted(handler.getTask(), trace) }
    }

    @Override
    void onProcessCached(TaskHandler handler, TraceRecord trace){
        log.trace "== event cached pid=${handler.getTask().processor.id}; status=$handler.status"
        agent.send { data.markCached(handler.getTask(), trace) }
    }

    WorkflowStats getStats() {
        return agent.getValue()
    }

    WorkflowStats getQuickStats() {
        agent.getQuickValue()
    }

    boolean hasProgressRecords() {
        return data.getProgressLength()
    }

    long getChangeTimestamp() {
        return data.changeTimestamp
    }
}
