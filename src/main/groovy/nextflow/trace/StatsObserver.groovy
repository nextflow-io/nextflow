package nextflow.trace

import nextflow.Session
import nextflow.processor.ErrorStrategy
import nextflow.processor.TaskHandler
import nextflow.processor.TaskProcessor
/**
 * Collect workflow execution statistics
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class StatsObserver implements TraceObserver {

    private WorkflowStats stats = new WorkflowStats()

    WorkflowStats getStats() { stats }


    protected long getCpuTime( TraceRecord record ) {
        // note: `realtime` field is not available when task metrics are enabled
        if( !record.containsKey('realtime') )
            return 0
        final time = (long)record.get('realtime')
        final cpus = (long)record.get('cpus')
        return time * cpus / 1000 as long
    }


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
        synchronized( stats ) {
            if( record.get('status') == 'COMPLETED' ) {
                stats.completed++
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
    }

    @Override
    void onProcessCached(TaskHandler handler) {
        final record = handler.getTraceRecord()
        synchronized (stats) {
            stats.timeCached += getCpuTime(record)
            stats.cached
        }
    }

    @Override
    boolean enableMetrics() {
        return false
    }
}
