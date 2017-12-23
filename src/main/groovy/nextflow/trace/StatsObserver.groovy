package nextflow.trace
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.processor.ErrorStrategy
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
        if( !record ) {
            log.debug "WARN: Unable to find trace record for task id=${handler.task?.id}"
            return
        }

        updateTasksCompleted(record)
    }

    @Override
    void onProcessCached(TaskHandler handler) {
        final record = handler.getTraceRecord()
        if( !record ) {
            log.debug "WARN: Unable to find trace record for task id=${handler.task?.id}"
            return
        }

        updateTasksCached(record)
    }

    @Override
    boolean enableMetrics() {
        return false
    }

    /**
     * Update succeed and failed task stats
     *
     * @param record The {@link TraceRecord} object representing the task metrics
     */
    protected synchronized void updateTasksCompleted(TraceRecord record) {
        if( record.get('status') == 'COMPLETED' ) {
            stats.succeed++
            stats.timeSucceed += getCpuTime(record)
        }
        else {
            stats.timeFailed += getCpuTime(record)

            final action = record.get('error_action')
            if( action == ErrorStrategy.IGNORE.toString() ) {
                stats.ignored++
            }
            else {
                stats.failed++
            }
        }
    }

    /**
     * Update cached task stats
     *
     * @param record The {@link TraceRecord} object representing the task metrics
     */
    protected synchronized void updateTasksCached(TraceRecord record) {
        stats.timeCached += getCpuTime(record)
        stats.cached++
    }

    /**
     * Computes the CPU-seconds used by the task
     *
     * @param record The {@link TraceRecord} object representing the task metrics
     * @return The CPU-seconds used by the specified task
     */
    protected long getCpuTime( TraceRecord record ) {
        // note: `realtime` field is not available when task metrics are not enabled
        if( !record.containsKey('realtime') )
            return 0
        final time = (long)record.get('realtime')
        final cpus = (long)record.get('cpus')
        return time * cpus / 1000 as long
    }

}
