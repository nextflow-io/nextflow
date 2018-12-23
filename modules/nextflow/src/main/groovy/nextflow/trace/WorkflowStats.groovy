/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.processor.ErrorStrategy
import nextflow.util.Duration
/**
 * Value object representing the workflow execution statistics
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@EqualsAndHashCode
@CompileStatic
class WorkflowStats {

    static private final DecimalFormat DECIMAL_FMT

    static private final DecimalFormat INTEGER_FMT

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

    private int succeedCount

    private int cachedCount

    private int failedCount

    private int ignoredCount

    String getSucceedCountFmt() {
        INTEGER_FMT.format(succeedCount)
    }

    String getCachedCountFmt() {
        INTEGER_FMT.format(cachedCount)
    }

    String getFailedCountFmt() {
        INTEGER_FMT.format(failedCount)
    }

    String getIgnoredCountFmt() {
        INTEGER_FMT.format(ignoredCount)
    }


    float getSucceedPct() {
        int tot = succeedCount + cachedCount + ignoredCount + failedCount
        tot ? Math.round(succeedCount / tot * 10000.0 as float) / 100.0 as float : 0
    }

    float getCachedPct() {
        def tot = succeedCount + cachedCount + ignoredCount + failedCount
        tot ? Math.round(cachedCount / tot * 10000.0 as float) / 100.0 as float : 0
    }

    float getIgnoredPct() {
        def tot = succeedCount + cachedCount + ignoredCount + failedCount
        tot ? Math.round(ignoredCount / tot * 10000.0 as float) / 100.0 as float : 0
    }

    float getFailedPct() {
        def tot = succeedCount + cachedCount + ignoredCount + failedCount
        tot ? Math.round(failedCount / tot * 10000.0 as float) / 100.0 as float : 0
    }

    /**
     * @return Overall workflow compute time (CPUs-seconds) for task executed successfully
     */
    Duration getSucceedDuration() { new Duration(succeedMillis) }

    /**
     * @return Overall workflow compute time (CPUs-seconds) for failed tasks
     */
    Duration getFailedDuration() { new Duration(failedMillis) }

    /**
     * @return Overall workflow compute time (CPUs-seconds) for cached tasks
     */
    Duration getCachedDuration() { new Duration(cachedMillis) }

    /**
     * @return Succeed tasks count
     */
    int getSucceedCount() { succeedCount }

    /**
     * @return Failed tasks count
     */
    int getFailedCount() { failedCount }

    /**
     * @return Ignored tasks count
     */
    int getIgnoredCount() { ignoredCount }

    /**
     * @return Cached tasks count
     */
    int getCachedCount() { cachedCount }

    /**
     * @return A formatted string representing the overall execution time as CPU-Hours
     */
    String getComputeTimeFmt() {

        final total = (succeedMillis + cachedMillis + failedMillis) / 1000
        if( total < 180 )
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

    /**
     * Update succeed and failed task stats
     *
     * @param record The {@link TraceRecord} object representing the task metrics
     */
    @PackageScope
    synchronized void updateTasksCompleted(TraceRecord record) {
        if( record.get('status') == 'COMPLETED' ) {
            succeedCount++
            succeedMillis += getCpuTime(record)
        }
        else {
            failedMillis += getCpuTime(record)

            final action = record.get('error_action')
            if( action == ErrorStrategy.IGNORE.toString() ) {
                ignoredCount++
            }
            else {
                failedCount++
            }
        }
    }

    /**
     * Update cached task stats
     *
     * @param record The {@link TraceRecord} object representing the task metrics
     */
    @PackageScope
    synchronized void updateTasksCached(TraceRecord record) {
        cachedMillis += getCpuTime(record)
        cachedCount++
    }

    String toString() {
        "WorkflowStats[" +
                "succeedCount=${getSucceedCount()}; " +
                "failedCount=${getFailedCount()}; " +
                "ignoredCount=${getIgnoredCount()}; " +
                "cachedCount=${getCachedCount()}; " +
                "succeedDuration=${getSucceedDuration()}; " +
                "failedDuration=${getFailedDuration()}; " +
                "cachedDuration=${getCachedDuration()}" +
                "]"
    }
}
