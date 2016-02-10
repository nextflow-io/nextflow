package nextflow.processor

import groovy.transform.CompileStatic

/**
 * Model a task execution fault
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class TaskFault {
    Throwable error
    String report
    TaskRun task
}
