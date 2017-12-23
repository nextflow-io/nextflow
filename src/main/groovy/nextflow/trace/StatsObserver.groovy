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
