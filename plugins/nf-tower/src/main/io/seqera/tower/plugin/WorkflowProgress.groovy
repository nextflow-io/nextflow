/*
 * Copyright (c) 2019, Seqera Labs.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * This Source Code Form is "Incompatible With Secondary Licenses", as
 * defined by the Mozilla Public License, v. 2.0.
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
