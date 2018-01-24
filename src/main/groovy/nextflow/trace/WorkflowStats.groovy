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

    private long succeedMillis

    private long cachedMillis

    private long failedMillis

    private int succeedCount

    private int cachedCount

    private int failedCount

    private int ignoredCount

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
    String getComputeTimeString() {

        final total = (succeedMillis + cachedMillis + failedMillis) / 1000
        if( total < 180 )
            return '(a few seconds)'

        final formatSymbols = new DecimalFormatSymbols();
        formatSymbols.setDecimalSeparator('.' as char);
        formatSymbols.setGroupingSeparator("'" as char);
        final fmt = new DecimalFormat("#,##0.0", formatSymbols)
        fmt.setRoundingMode(RoundingMode.HALF_UP)

        def result = fmt.format(total/3600)
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
