package nextflow.processor

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
interface TaskMonitor {

    /**
     * Put a {@code TaskHandler} instance into the queue of tasks to be processed.
     * Subclasses provide concrete logic to manage them.
     *
     * Note: depending the task scheduling implementation this operation may be blocking
     * to await there's free space in the tasks queue
     *
     * @param handler {@code TaskHandler} instance
     */
    def void put(TaskHandler handler)

    /**
     * Remove the {@code TaskHandler} instance from the queue of tasks ro be processed
     *
     * @param handler A not null {@code TaskHandler} instance
     */
    def boolean remove(TaskHandler handler)

    /**
     * Start the monitoring activity for the queued tasks
     * @return The instance itself, useful to chain methods invocation
     */
    def TaskMonitor start()
}
