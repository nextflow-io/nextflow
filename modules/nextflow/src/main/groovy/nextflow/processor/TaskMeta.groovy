package nextflow.processor

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.transform.ToString

/**
 * Model task execution meta used for testflow tests
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@ToString(includeNames = true, includePackage = false)
class TaskMeta {
    final int index
    final String name
    final String hash
    final Path workDir
    final String tag

    TaskMeta(TaskRun task) {
        index = task.index as int
        name = task.name
        hash = task.hash.toString()
        workDir = task.workDir
        tag = task.config.tag
    }
}
