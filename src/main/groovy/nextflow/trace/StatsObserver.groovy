/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
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
    void onProcessSubmit(TaskHandler handler) {

    }

    @Override
    void onProcessStart(TaskHandler handler) {

    }

    @Override
    void onProcessComplete(TaskHandler handler) {
        final record = handler.getTraceRecord()
        if( record ) {
            stats.updateTasksCompleted(record)
        }
        else {
            log.debug "WARN: Unable to find trace record for task id=${handler.task?.id}"
        }
    }

    @Override
    void onProcessCached(TaskHandler handler) {
        final record = handler.getTraceRecord()
        if( record ) {
            stats.updateTasksCached(record)
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
