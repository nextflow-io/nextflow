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
/**
 * Collect workflow execution progress and statistics
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class StatsObserver implements TraceObserver, ProgressState {

    private WorkflowStats stats = new WorkflowStats()

    private Session session

    private List<ProgressRecord> records = new ArrayList<>(100)

    private volatile long changeTimestamp

    WorkflowStats getStats() { stats }

    StatsObserver(Session session) {
        this.session = session
    }

    @Override
    void onFlowComplete() {
        log.debug "Workflow completed > $stats"
    }

    @Override
    void onProcessCreate(TaskProcessor process){
        records[ process.id-1 ] = new ProgressRecord(process.id, process.name)
        changeTimestamp = System.currentTimeMillis()
    }

    @Override
    void onProcessTerminate( TaskProcessor processor ) {
        final state = records[ processor.id-1 ]
        state?.markTerminated()
        changeTimestamp = System.currentTimeMillis()
    }

    @Override
    void onProcessPending(TaskHandler handler, TraceRecord trace){
        final state = records[ handler.getTask().processor.id-1 ]
        state?.markPending()
        changeTimestamp = System.currentTimeMillis()
    }

    @Override
    void onProcessSubmit(TaskHandler handler, TraceRecord trace){
        final state = records[ handler.getTask().processor.id-1 ]
        state?.markSubmitted(handler.getTask())
        changeTimestamp = System.currentTimeMillis()
    }

    @Override
    void onProcessStart(TaskHandler handler, TraceRecord trace){
        final state = records[ handler.getTask().processor.id-1 ]
        state?.markRunning(handler.getTask())
        changeTimestamp = System.currentTimeMillis()
    }

    @Override
    void onProcessComplete(TaskHandler handler, TraceRecord trace) {
        ProgressRecord state = records[ handler.getTask().processor.id-1 ]
        state?.markComplete(handler.getTask())
        changeTimestamp = System.currentTimeMillis()

        // update stats
        if( trace ) {
            stats.updateTasksCompleted(trace)
        }
        else {
            log.debug "WARN: Unable to find trace record for task id=${handler.task?.id}"
        }
    }

    @Override
    void onProcessCached(TaskHandler handler, TraceRecord trace){
        final state = records[ handler.getTask().processor.id-1 ]
        state?.markCached(handler.getTask(), trace)
        changeTimestamp = System.currentTimeMillis()
        // update stats
        if( trace ) {
            stats.updateTasksCached(trace)
        }
    }

    @Override
    List<ProgressRecord> getProgress() {
        def result = new ArrayList<ProgressRecord>(records.size())
        for( int i=0; i<records.size(); i++ ) {
            result.add(i, records.get(i).clone())
        }
        return result
    }

    @Override
    int getProgressLength() {
        records.size()
    }

    @Override
    long getChangeTimestamp() {
        changeTimestamp
    }

}
