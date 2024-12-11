package nextflow.executor

import groovy.transform.CompileStatic
import nextflow.processor.TaskHandler
import nextflow.processor.TaskRun
import nextflow.trace.TraceRecord

/**
 * Implements a {@link TaskHandler} instance for nextflow stored task ie.
 * tasks whose execution is skipped due the use of the `storeDir` directive.
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class StoredTaskHandler extends TaskHandler {

    StoredTaskHandler(TaskRun task) {
        super(task)
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
    protected void killTask() {
        throw new UnsupportedOperationException()
    }

    @Override
    void submit() {
        throw new UnsupportedOperationException()
    }

    @Override
    String getStatusString() {
        "STORED"
    }

    /**
     * @return Stored tasks are not supposed to have a trace record, therefore returns {@code null}
     */
    @Override
    TraceRecord getTraceRecord() {
        return null
    }
}
