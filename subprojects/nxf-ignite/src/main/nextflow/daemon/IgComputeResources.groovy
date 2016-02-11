package nextflow.daemon

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString
import nextflow.processor.TaskRun
import nextflow.util.Duration
import nextflow.util.MemoryUnit

/**
 * Model the computing resources required by a task
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@ToString
@EqualsAndHashCode
@CompileStatic
class IgComputeResources implements Serializable {

    int cpus

    MemoryUnit memory

    MemoryUnit disk

    Duration time

    IgComputeResources() {}

    IgComputeResources( TaskRun task ) {
        cpus = task.config.getCpus()
        memory = task.config.getMemory()
        disk = task.config.getDisk()
        time = task.config.getTime()
    }

}
