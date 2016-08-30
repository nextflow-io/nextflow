package nextflow.executor

import nextflow.processor.TaskHandler
import nextflow.processor.TaskRun
import nextflow.trace.TraceRecord

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CachedTaskHandler extends TaskHandler {

    private TraceRecord trace

    CachedTaskHandler(TaskRun task, TraceRecord trace) {
        super(task)
        this.trace = trace
        this.trace.setCached(true)
    }

    @Override
    boolean checkIfRunning() {
        return false
    }

    @Override
    boolean checkIfCompleted() {
        return true
    }

    @Override
    void kill() {
        throw new UnsupportedOperationException()
    }

    @Override
    void submit() {
        throw new UnsupportedOperationException()
    }

    @Override
    String getStatusString() {
        "CACHED"
    }

    @Override
    TraceRecord getTraceRecord() {
        return trace
    }

}
