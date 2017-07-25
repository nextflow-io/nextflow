package nextflow.executor

import java.nio.file.Path

import groovy.transform.CompileStatic
import io.swagger.client.api.TaskServiceApi
import io.swagger.client.model.Body
import io.swagger.client.model.OUTPUTONLYExecutors
import nextflow.processor.TaskHandler
import nextflow.processor.TaskMonitor
import nextflow.processor.TaskPollingMonitor
import nextflow.processor.TaskRun
import nextflow.util.Duration
import nextflow.util.Escape

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class TesExecutor extends Executor {


    /**
     * Create a a queue holder for this executor
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

        // create the wrapper script
        final bash = new BashWrapperBuilder(task)
        bash.build()

        new TesTaskHandler(task, this)
    }
}


@CompileStatic
class TesTaskHandler extends TaskHandler {

    final TesExecutor executor

    private final Path exitFile

    private final Long wallTimeMillis

    private final Path wrapperFile

    private final Path outputFile

    private final Path errorFile

    private TaskServiceApi api

    private String requestId

    TesTaskHandler(TaskRun task, TesExecutor executor) {
        super(task)
        this.executor = executor
        this.api = new TaskServiceApi()
        this.exitFile = task.workDir.resolve(TaskRun.CMD_EXIT)
        this.outputFile = task.workDir.resolve(TaskRun.CMD_OUTFILE)
        this.errorFile = task.workDir.resolve(TaskRun.CMD_ERRFILE)
        this.wrapperFile = task.workDir.resolve(TaskRun.CMD_RUN)
    }

    @Override
    boolean checkIfRunning() {
        return false
    }

    @Override
    boolean checkIfCompleted() {
        return false
    }

    @Override
    void kill() {

    }

    @Override
    void submit() {

        // the cmd list to launch it
        def job = new ArrayList(BashWrapperBuilder.BASH) << wrapperFile.getName()
        List cmd = ['/bin/bash','-c', job.join(' ') + " &> ${Escape.path(task.workDir)}/${TaskRun.CMD_LOG}" ]

        def exec = new OUTPUTONLYExecutors()
        exec.setCmd(cmd)

        final body = new Body()
        body.executors = [exec]
        body.inputs = []
        body.outputs = []

        final task = api.createTask(body)
        requestId = task.id
    }
}
