package nextflow.trace
import java.text.DecimalFormat

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

        def fmt = new DecimalFormat("0.#")
        def total = (succeedMillis + cachedMillis + failedMillis) / 1000
        if( total < 180 )
            return '(a few seconds)'

        def result = String.format('%.1f', total/3600)
        if( cachedMillis || failedMillis ) {
            result += ' ('
            def items = []
            if( cachedMillis ) {
                items << fmt.format(cachedMillis/10/total) + '% cached'
            }
            if( failedMillis ) {
                items << fmt.format(failedMillis/10/total) + '% failed'
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
        // note: `realtime` field is not available when task metrics are not enabled
        if( !record.containsKey('realtime') )
            return 0
        final time = (long)record.get('realtime')
        final cpus = (long)record.get('cpus')
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
