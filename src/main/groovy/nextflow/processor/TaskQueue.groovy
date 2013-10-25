package nextflow.processor

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
interface TaskQueue {

    def void addTask(TaskHandler handler)


}
