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

import java.math.RoundingMode
import java.text.DecimalFormat
import java.text.DecimalFormatSymbols

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.util.logging.Slf4j
import nextflow.processor.ErrorStrategy
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import nextflow.util.Duration
import nextflow.util.MemoryUnit
/**
 * Value object representing the workflow execution statistics
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@EqualsAndHashCode
@CompileStatic
class WorkflowStats implements Cloneable {

    static public final short MIN_SECS = 180

    static private final DecimalFormat DECIMAL_FMT

    static private final DecimalFormat INTEGER_FMT

    private Map<Integer, ProgressRecord> records = new TreeMap<>()

    private transient volatile long changeTimestamp

    static {
        final formatSymbols = new DecimalFormatSymbols()
        formatSymbols.setDecimalSeparator('.' as char)
        formatSymbols.setGroupingSeparator("'" as char)

        DECIMAL_FMT = new DecimalFormat("#,##0.0", formatSymbols)
        DECIMAL_FMT.setRoundingMode(RoundingMode.HALF_UP)

        INTEGER_FMT = new DecimalFormat("#,##0", formatSymbols)
        INTEGER_FMT.setRoundingMode(RoundingMode.HALF_UP)
    }

    private long succeedMillis
    private long cachedMillis
    private long failedMillis

    private int succeededCount
    private int cachedCount
    private int failedCount
    private int ignoredCount
    private int pendingCount
    private int submittedCount
    private int runningCount
    private int retriesCount
    private int abortedCount

    private int loadCpus
    private long loadMemory
    private int peakRunning
    private int peakCpus
    private long peakMemory

    WorkflowStats clone() {
        final result = (WorkflowStats)super.clone()
        result.records = new TreeMap<>(this.records)
        for( Integer id : records.keySet() ) {
            result.records[id] = records[id].clone()
        }
        return result
    }

    static private long gtz(long x) {
        x >= 0 ? x : 0
    }

    static private int gtz(int x) {
        x >= 0 ? x : 0
    }

    String getSucceedCountFmt() {
        INTEGER_FMT.format(gtz(succeededCount))
    }

    String getCachedCountFmt() {
        INTEGER_FMT.format(gtz(cachedCount))
    }

    String getFailedCountFmt() {
        INTEGER_FMT.format(gtz(failedCount))
    }

    String getIgnoredCountFmt() {
        INTEGER_FMT.format(gtz(ignoredCount))
    }

    float getSucceedPct() {
        int tot = gtz(succeededCount + cachedCount + ignoredCount + failedCount)
        tot ? Math.round(succeededCount / tot * 10000.0 as float) / 100.0 as float : 0
    }

    float getCachedPct() {
        def tot = gtz(succeededCount + cachedCount + ignoredCount + failedCount)
        tot ? Math.round(gtz(cachedCount) / tot * 10000.0 as float) / 100.0 as float : 0
    }

    float getIgnoredPct() {
        def tot = gtz(succeededCount + cachedCount + ignoredCount + failedCount)
        tot ? Math.round(gtz(ignoredCount) / tot * 10000.0 as float) / 100.0 as float : 0
    }

    float getFailedPct() {
        def tot = gtz(succeededCount + cachedCount + ignoredCount + failedCount)
        tot ? Math.round(gtz(failedCount) / tot * 10000.0 as float) / 100.0 as float : 0
    }

    protected Duration makeDuration(long value) {
        new Duration(gtz(value))
    }

    /**
     * @return Overall workflow compute time (CPUs-seconds) for task executed successfully
     */
    Duration getSucceedDuration() { makeDuration(succeedMillis) }

    /**
     * @return Overall workflow compute time (CPUs-seconds) for failed tasks
     */
    Duration getFailedDuration() { makeDuration(failedMillis) }

    /**
     * @return Overall workflow compute time (CPUs-seconds) for cached tasks
     */
    Duration getCachedDuration() { makeDuration(cachedMillis) }

    /**
     * @return Succeed tasks count
     */
    @Deprecated
    int getSucceedCount() { gtz(succeededCount) }

    /**
     * @return Succeed tasks count
     */
    int getSucceededCount() { gtz(succeededCount) }

    /**
     * @return Failed tasks count
     */
    int getFailedCount() { gtz(failedCount) }

    /**
     * @return Ignored tasks count
     */
    int getIgnoredCount() { gtz(ignoredCount) }

    /**
     * @return Cached tasks count
     */
    int getCachedCount() { gtz(cachedCount) }

    int getPendingCount() { gtz(pendingCount) }

    int getSubmittedCount() { gtz(submittedCount) }

    int getRunningCount() { gtz(runningCount) }

    int getRetriesCount() { gtz(retriesCount) }

    int getAbortedCount() { gtz(abortedCount) }

    int getLoadCpus() { gtz(loadCpus) }

    long getLoadMemory() { gtz(loadMemory) }

    String getLoadMemoryFmt() { MemoryUnit.of(getLoadMemory()).toString() }

    int getPeakRunning() { gtz(peakRunning) }

    long getPeakCpus() { gtz(peakCpus) }

    long getPeakMemory() { gtz(peakMemory) }

    String getPeakMemoryFmt() { MemoryUnit.of(getPeakMemory()).toString() }

    /**
     * @return A formatted string representing the overall execution time as CPU-Hours
     */
    String getComputeTimeFmt() {

        final total = gtz(succeedMillis + cachedMillis + failedMillis) / 1000
        if( total < MIN_SECS )
            return '(a few seconds)'


        def result = DECIMAL_FMT.format(total/3600)
        if( cachedMillis || failedMillis ) {
            final perc = new DecimalFormat("0.#")
            result += ' ('
            def items = []
            if( cachedMillis ) {
                items << perc.format(cachedMillis/10/total) + '% cached'
            }
            if( failedMillis ) {
                items << perc.format(failedMillis/10/total) + '% failed'
            }
            result += items.join(', ')
            result += ')'
        }

        return result
    }


    /**
     * Computes the CPU-millis used by the task
     *
     * @param record The {@link TraceRecord} object representing the task metrics
     * @return The CPU-millis used by the specified task
     */
    protected long getCpuTime( TraceRecord record ) {
        final time = (long) (record.get('realtime') ?: 0)
        final cpus = (int) (record.get('cpus') ?: 1)
        return time * cpus
    }

    String toString() {
        "WorkflowStats[" +
                "succeededCount=${getSucceededCount()}; " +
                "failedCount=${getFailedCount()}; " +
                "ignoredCount=${getIgnoredCount()}; " +
                "cachedCount=${getCachedCount()}; " +
                "pendingCount=${pendingCount}; " +
                "submittedCount=${submittedCount}; " +
                "runningCount=${runningCount}; " +
                "retriesCount=${retriesCount}; " +
                "abortedCount=${abortedCount}; " +
                "succeedDuration=${getSucceedDuration()}; " +
                "failedDuration=${getFailedDuration()}; " +
                "cachedDuration=${getCachedDuration()};" +
                "loadCpus=${loadCpus}; " +
                "loadMemory=${getLoadMemoryFmt()}; " +
                "peakRunning=${peakRunning}; " +
                "peakCpus=${peakCpus}; " +
                "peakMemory=${getPeakMemoryFmt()}; " +
                "]"
    }

    private ProgressRecord getOrCreateRecord(TaskProcessor process) {
        def pid = process.getId()
        def rec = records.get(pid)
        if( !rec ) {
            rec = new ProgressRecord(pid, process.getName())
            records.put(pid, rec)
        }
        return rec
    }

    void markCreated(TaskProcessor process) {
        getOrCreateRecord(process)
        changeTimestamp = System.currentTimeMillis()
    }

    void markPending(TaskProcessor process) {
        final state = getOrCreateRecord(process)
        state.pending ++
        // global counters
        this.pendingCount ++
        this.changeTimestamp = System.currentTimeMillis()
    }

    void markSubmitted(TaskRun task) {
        final state = getOrCreateRecord(task.processor)
        state.hash = task.hashLog
        state.taskName = task.name
        state.pending --
        state.submitted ++
        // global counters
        this.pendingCount --
        this.submittedCount ++
        this.changeTimestamp = System.currentTimeMillis()
    }

    void markRunning(TaskRun task) {
        final state = getOrCreateRecord(task.processor)
        state.with {
            submitted --
            running ++
            // update current load
            loadCpus += task.getConfig().getCpus()
            loadMemory += (task.getConfig().getMemory()?.toBytes() ?: 0)
            // update peaks
            if( peakRunning < running )
                peakRunning = running
            if( peakCpus < loadCpus )
                peakCpus = loadCpus
            if( peakMemory < loadMemory )
                peakMemory = loadMemory
        }

        // global counters
        submittedCount --
        runningCount ++
        // update current load
        loadCpus += task.getConfig().getCpus()
        loadMemory += (task.getConfig().getMemory()?.toBytes() ?: 0)

        // update peaks
        if( peakRunning < runningCount )
            peakRunning = runningCount
        if( peakCpus < loadCpus )
            peakCpus = loadCpus
        if( peakMemory < loadMemory )
            peakMemory = loadMemory

        changeTimestamp = System.currentTimeMillis()

    }

    void markCompleted(TaskRun task, TraceRecord trace) {
        ProgressRecord state = getOrCreateRecord(task.processor)
        state.taskName = task.name
        state.hash = task.hashLog
        state.running --
        state.loadCpus -= task.getConfig().getCpus()
        state.loadMemory -= (task.getConfig().getMemory()?.toBytes() ?: 0)

        this.runningCount --
        this.loadCpus -= task.getConfig().getCpus()
        this.loadMemory -= (task.getConfig().getMemory()?.toBytes() ?: 0)

        if( task.failed ) {
            state.failed ++
            this.failedCount ++
            this.failedMillis += trace ? getCpuTime(trace) : 0

            if( task.errorAction == ErrorStrategy.RETRY ) {
                state.retries ++
                this.retriesCount ++
            }

            else if( task.errorAction == ErrorStrategy.IGNORE ) {
                state.ignored ++
                this.ignoredCount ++
            }

            else if( !task.errorAction?.soft )
                state.errored |= true
        }
        else if( task.aborted ) {
            state.aborted ++
            this.abortedCount ++
        }
        else {
            state.succeeded ++
            this.succeededCount++
            this.succeedMillis += trace ? getCpuTime(trace) : 0
        }

        changeTimestamp = System.currentTimeMillis()
    }

    void markCached(TaskRun task, TraceRecord trace) {
        final state = getOrCreateRecord(task.processor)
        if( trace ) {
            state.cached++
            state.hash = task.hashLog
            state.taskName = task.name
            // global counters
            this.cachedMillis += getCpuTime(trace)
            this.cachedCount++
        }
        else {
            state.stored++
            state.hash = 'skipped'
            state.taskName = task.name
        }
        changeTimestamp = System.currentTimeMillis()
    }

    void markTerminated(TaskProcessor processor) {
        final state = getOrCreateRecord(processor)
        state.terminated = true
        changeTimestamp = System.currentTimeMillis()
    }

    List<ProgressRecord> getProcesses() {
        new ArrayList<ProgressRecord>(records.values())
    }

    int getProgressLength() {
        records.size()
    }

    long getChangeTimestamp() {
        changeTimestamp
    }

}
