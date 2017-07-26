package nextflow.executor
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.processor.TaskHandler
import nextflow.processor.TaskMonitor
import nextflow.processor.TaskPollingMonitor
import nextflow.processor.TaskRun
import nextflow.util.Duration
import nextflow.util.ServiceName
/**
 * Experimental TES executor
 *
 * See https://github.com/ga4gh/task-execution-schemas/
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
@ServiceName('tes')
class TesExecutor extends Executor {

    /**
     * Create a a queue holder for this executor
     *
     * @return
     */
    def TaskMonitor createTaskMonitor() {
        return TaskPollingMonitor.create(session, name, 50, Duration.of('1 sec'))
    }


    /*
     * Prepare and launch the task in the underlying execution platform
     */
    @Override
    TaskHandler createTaskHandler(TaskRun task) {
        assert task
        assert task.workDir

        log.debug "Launching process > ${task.name} -- work folder: ${task.workDir}"

        new TesTaskHandler(task, this)
    }
}



