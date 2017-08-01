package nextflow.executor

import static nextflow.processor.TaskStatus.COMPLETED
import static nextflow.processor.TaskStatus.RUNNING

import java.nio.file.Path
import java.nio.file.Paths

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.ga4gh.tes.client.api.TaskServiceApi
import nextflow.ga4gh.tes.client.model.Body
import nextflow.ga4gh.tes.client.model.OUTPUTONLYExecutors
import nextflow.ga4gh.tes.client.model.OUTPUTONLYInputs
import nextflow.ga4gh.tes.client.model.OUTPUTONLYTasks
import nextflow.processor.TaskHandler
import nextflow.processor.TaskRun
import nextflow.processor.TaskStatus

@Slf4j
@CompileStatic
class TesTaskHandler extends TaskHandler {

    static final Path containerWorkDir = Paths.get('/work')

    final List<OUTPUTONLYTasks.StateEnum> COMPLETE_STATUSES = [OUTPUTONLYTasks.StateEnum.COMPLETE, OUTPUTONLYTasks.StateEnum.ERROR, OUTPUTONLYTasks.StateEnum.SYSTEM_ERROR, OUTPUTONLYTasks.StateEnum.CANCELED]

    final List<OUTPUTONLYTasks.StateEnum> STARTED_STATUSES = [OUTPUTONLYTasks.StateEnum.RUNNING, OUTPUTONLYTasks.StateEnum.PAUSED] + COMPLETE_STATUSES

    final TesExecutor executor

    private final Path exitFile

    private final Path wrapperFile

    private final Path outputFile

    private final Path errorFile

    private final Path logFile

    private final Path envFile

    private final Path scriptFile

    private final Path inputFile

    private final Path stubFile

    private final Path traceFile

    private TaskServiceApi api

    private String requestId

    TesTaskHandler(TaskRun task, TesExecutor executor) {
        super(task)
        this.executor = executor
        this.api = new TaskServiceApi()

        this.logFile = task.workDir.resolve(TaskRun.CMD_LOG)
        this.envFile = task.workDir.resolve(TaskRun.CMD_ENV)
        this.scriptFile = task.workDir.resolve(TaskRun.CMD_SCRIPT)
        this.inputFile =  task.workDir.resolve(TaskRun.CMD_INFILE)
        this.outputFile = task.workDir.resolve(TaskRun.CMD_OUTFILE)
        this.errorFile = task.workDir.resolve(TaskRun.CMD_ERRFILE)
        this.exitFile = task.workDir.resolve(TaskRun.CMD_EXIT)
        this.wrapperFile = task.workDir.resolve(TaskRun.CMD_RUN)
        this.stubFile = task.workDir.resolve(TaskRun.CMD_STUB)
        this.traceFile = task.workDir.resolve(TaskRun.CMD_TRACE)

    }

    @Override
    boolean checkIfRunning() {

        if( requestId && isSubmitted() ) {
            final response = api.getTask(requestId, null)
            final started = response.state in STARTED_STATUSES
            if( started ) {
                log.debug "Task started > $task.name"
                status = RUNNING
                return true
            }
        }

        return false
    }

    @Override
    boolean checkIfCompleted() {
        if( !isRunning() ) {
            return false
        }

        final response = api.getTask(requestId, null)
        if( response.state in COMPLETE_STATUSES ) {
            // finalize the task
            log.debug "Task completed > $task.name"
            task.exitStatus = readExitFile()
            task.stdout = outputFile
            task.stderr = errorFile
            status = COMPLETED
            return true
        }

        return false
    }

    private int readExitFile() {
        try {
            exitFile.text as Integer
        }
        catch( Exception e ) {
            log.debug "Cannot read exitstatus for task: `$task.name`", e
            return Integer.MAX_VALUE
        }
    }

    @Override
    void kill() {
        if( requestId )
            api.cancelTask(requestId)
        else
            log.debug "Invalid kill request -- missing requestId"
    }

    @Override
    void submit() {

        // create task wrapper
        final bash = new TesBashBuilder(task)
        bash.build()

        // the cmd list to launch it
        def job = new ArrayList(BashWrapperBuilder.BASH) << wrapperFile.getName()
        List cmd = ['/bin/bash','-c', job.join(' ') + " &> $TaskRun.CMD_LOG" ]

        def exec = new OUTPUTONLYExecutors()
        exec.cmd = cmd
        exec.imageName = task.container
        exec.workdir = containerWorkDir

        final body = new Body()

        // add task control files
        body.addInputsItem(inItem(scriptFile))
        body.addInputsItem(inItem(wrapperFile))

        // add task input files
        if(envFile.exists()) body.addInputsItem(inItem(envFile))
        if(inputFile.exists()) body.addInputsItem(inItem(inputFile))
        if(stubFile.exists()) body.addInputsItem(inItem(stubFile))

        task.getInputFilesMap().each { String name, Path path ->
            body.addInputsItem(inItem(path,name))
        }

        // add the task output files
        body.addOutputsItem(outItem(outputFile.name))
        body.addOutputsItem(outItem(errorFile.name))
        body.addOutputsItem(outItem(logFile.name))
        body.addOutputsItem(outItem(exitFile.name))

        task.outputFilesNames.each { fileName ->
            body.addOutputsItem(outItem(fileName))
        }

        // add the executor
        body.executors = [exec]

        // submit the task
        final task = api.createTask(body)
        requestId = task.id
        status = TaskStatus.SUBMITTED
    }


    private OUTPUTONLYInputs inItem( Path realPath, String fileName = null) {
        def result = new OUTPUTONLYInputs()
        result.url = realPath.toUri().toString()
        result.path = fileName ? containerWorkDir.resolve(fileName) : containerWorkDir.resolve(realPath.getName())
        log.debug "Adding INPUT file: $result"
        return result
    }

    private OUTPUTONLYInputs outItem( String fileName ) {
        def result = new OUTPUTONLYInputs()
        result.path = containerWorkDir.resolve(fileName)
        result.url = task.workDir.resolve(fileName).uri.toString()
        log.debug "Adding OUTPUT file: $result"
        return result
    }

}