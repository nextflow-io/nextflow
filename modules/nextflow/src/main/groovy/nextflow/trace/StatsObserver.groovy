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
 * Collect workflow execution statistics
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class StatsObserver implements TraceObserver {

    private WorkflowStats stats = new WorkflowStats()

    WorkflowStats getStats() { stats }

    @Override
    void onFlowStart(Session session) {

    }

    @Override
    void onFlowComplete() {
        log.debug "Workflow completed > $stats"
    }

    @Override
    void onProcessCreate(TaskProcessor process) {

    }

    @Override
    void onProcessSubmit(TaskHandler handler, TraceRecord trace) {

    }

    @Override
    void onProcessStart(TaskHandler handler, TraceRecord trace) {

    }

    @Override
    void onProcessComplete(TaskHandler handler, TraceRecord trace) {
        if( trace ) {
            stats.updateTasksCompleted(trace)
        }
        else {
            log.debug "WARN: Unable to find trace record for task id=${handler.task?.id}"
        }
    }

    @Override
    void onProcessCached(TaskHandler handler, TraceRecord trace) {
        if( trace ) {
            stats.updateTasksCached(trace)
        }
        else {
            log.debug "WARN: Unable to find trace record for task id=${handler.task?.id}"
        }
    }

    @Override
    boolean enableMetrics() {
        return false
    }

}
