package nextflow.processor

import groovy.transform.Canonical
import nextflow.trace.TraceRecord

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Canonical
class TaskEntry {

    TraceRecord trace

    TaskContext context

}
