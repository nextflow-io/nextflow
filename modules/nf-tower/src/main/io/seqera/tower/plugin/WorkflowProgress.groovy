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

package io.seqera.tower.plugin

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import nextflow.trace.ProgressRecord
import nextflow.trace.WorkflowStats
/**
 * Simple facade to format workflow progress payload
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@EqualsAndHashCode
@CompileStatic
class WorkflowProgress {
    private WorkflowStats stats

    WorkflowProgress(WorkflowStats stats) {
        this.stats = stats
    }

    int getSucceeded() { stats.succeededCount }

    int getFailed() { stats.failedCount }

    int getIgnored() { stats.ignoredCount }

    int getCached() { stats.cachedCount }

    int getPending() { stats.pendingCount }

    int getSubmitted() { stats.submittedCount }

    int getRunning() { stats.runningCount }

    int getRetries() { stats.retriesCount }

    int getAborted() { stats.abortedCount }

    int getLoadCpus() { stats.loadCpus }

    long getLoadMemory() { stats.loadMemory }

    int getPeakRunning() { stats.peakRunning }

    long getPeakCpus() { stats.peakCpus }

    long getPeakMemory() { stats.peakMemory }

    List<ProgressRecord> getProcesses() { stats.getProcesses() }

}
