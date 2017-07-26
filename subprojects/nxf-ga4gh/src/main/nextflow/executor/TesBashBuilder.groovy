package nextflow.executor

import groovy.transform.CompileStatic
import nextflow.processor.TaskBean
import nextflow.processor.TaskRun

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class TesBashBuilder extends BashWrapperBuilder {

    TesBashBuilder(TaskRun task) {
        super(new TaskBean(task), new NoopFileCopyStrategy())
    }

    TesBashBuilder(TaskBean task) {
        super(task, new NoopFileCopyStrategy())
    }

    protected boolean containerInit() { false }
}
