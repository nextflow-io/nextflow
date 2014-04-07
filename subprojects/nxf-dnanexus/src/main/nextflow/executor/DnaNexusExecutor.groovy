/*
 * Copyright (c) 2013-2014, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2014, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
 */

package nextflow.executor
import static nextflow.processor.TaskHandler.Status.COMPLETED
import static nextflow.processor.TaskHandler.Status.RUNNING
import static nextflow.processor.TaskHandler.Status.SUBMITTED

import java.nio.file.Files
import java.nio.file.Path

import groovy.util.logging.Slf4j
import nextflow.fs.dx.DxFileSystem
import nextflow.fs.dx.DxFileSystemSerializer
import nextflow.fs.dx.DxPath
import nextflow.fs.dx.DxPathSerializer
import nextflow.fs.dx.api.DxApi
import nextflow.processor.TaskConfig
import nextflow.processor.TaskHandler
import nextflow.processor.TaskMonitor
import nextflow.processor.TaskPollingMonitor
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import nextflow.util.Duration
import nextflow.util.KryoHelper

/**
 * Executes script.nf indicated in dxapp.sh in the DnaNexus environment
 *
 * See https://www.dnanexus.com/
 *
 * @author Beatriz Martin San Juan <bmsanjuan@gmail.com>
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */


@Slf4j
@ServiceName('dnanexus')
class DnaNexusExecutor extends Executor {

    def void init() {
        super.init()
        // -- register the serializer for the Dx file system objects
        KryoHelper.register(DxPath, DxPathSerializer)
        KryoHelper.register(DxFileSystem, DxFileSystemSerializer)
    }

    def TaskMonitor createTaskMonitor() {
        return TaskPollingMonitor.create(session, name, 10, Duration.of('15 sec'))
    }

    /**
     * Returns the output of the task.
     * @param task
     * @return task.output
     */
    @Override
    TaskHandler createTaskHandler(TaskRun task) {

        /*
         * Setting the work directory
         */
        final scratch = task.workDirectory
        log.debug "Lauching process > ${task.name} -- work folder: $scratch"

        /*
         * Saving the environment to a file.
         */
        def taskEnvFile = null
        Map environment = task.processor.getProcessEnvironment()
        environment.putAll( task.getInputEnvironment() )
        final envBashText = TaskProcessor.bashEnvironmentScript(environment)
        // make sure to upload the file only there is an environment content
        if( envBashText ) {
            taskEnvFile = task.getCmdEnvironmentFile()
            taskEnvFile.text = envBashText
        }

        /*
         * In case there's a task input file, save it file
         */
        Path taskInputFile = null
        if (task.stdin) {
            taskInputFile = task.getCmdInputFile()
            taskInputFile.text = task.stdin
        }

        /*
         * Saving the task script file in the appropriate task's working folder
         */
        Path taskScriptFile = task.getCmdScriptFile()
        taskScriptFile.text = task.processor.normalizeScript(task.script.toString())
        log.debug "Creating script file for process > ${task.name}\n\n "

        /*
         * create the stage inputs script
         */
        def inputFiles = task.getInputFiles()
        def stageScript = stagingFilesScript( inputFiles, '; ' )

        /*
         * the DnaNexus job params object
         */
        def obj = [:]
        obj.task_name = task.name
        obj.task_script = (taskScriptFile as DxPath).getFileId()
        if( taskEnvFile ) {
            obj.task_env = (taskEnvFile as DxPath).getFileId()
        }
        if( taskInputFile ) {
            obj.task_input = (taskInputFile as DxPath).getFileId()
        }
        obj.stage_inputs = stageScript
        obj.output_files = task.getOutputFilesNames()

        def target = task.getTargetDir().toString()
        if( target.startsWith('dxfs://') )
            target = target.substring('dxfs://'.size())
        obj.target_dir = target

        new DxTaskHandler(task, taskConfig, this, obj)
    }


    @Override
    String stageInputFileScript( Path path, String target ) {
        if( !(path instanceof DxPath) ) {
            return super.stageInputFileScript(path,target)
        }

        if( Files.isDirectory(path) ) {
            def origin = path.toAbsolutePath().normalize();
            return "tmp=\$(mktemp -d); mkdir -p '$target'; dx download --no-progress -r ${origin} -o \$tmp; mv \$tmp/${origin}/* '${target}'"
        }
        else {
            return "dx download --no-progress ${(path as DxPath).getFileId()} -o $target"
        }

    }


    /**
     * Building the ObjectNode which will be set in the job.
     * Depending on whether we have the already checked task's input file or not,
     *
     * @param inputObj
     *          List formed by all the names and ids of the inputs declared.
     * @param outputs
     *          List formed by all the names of the outputs declared.
     * @param taskName
     *          String with the name of the task.
     * @param env
     *          String created with all the variables from it
     * @param scriptId
     *          Id of the task's script file
     * @param taskInputId
     *          (Compulsory if we have the named file) Id of the task's input file if we have a task's input; null if not.
     * @param instanceType
     *          Value of the instance type to be used (optional)
     */

    def Map createInputObject( Map inputObj, String instanceType ) {
        log.debug "Creating inputObj with: ${inputObj}"

        def root = [:]

        if(instanceType){
            def process = [ instanceType: instanceType ]
            root.systemRequirements = [process: process]
        }

        root.input = inputObj
        root.function = "process"

        return root

    }

}

/**
 * Handle a job execution in the DnaNexus platform
 */
@Slf4j
class DxTaskHandler extends TaskHandler {

    final Map inputParams

    final DnaNexusExecutor executor

    final private DxApi api

    private String processJobId

    private long lastStatusMillis

    private Map lastStatusResult



    protected DxTaskHandler(TaskRun task, TaskConfig config, DnaNexusExecutor executor, Map params, DxApi api = null ) {
        super(task, config)
        this.inputParams = params
        this.executor = executor
        this.api = api ?: DxApi.instance
    }

    @Override
    void submit() {

        // create the input parameters for the job to be executed
        def processJobInputHash = executor.createInputObject( inputParams, (String)taskConfig.instanceType )
        log.debug "New job parameters: ${processJobInputHash}"

        // Launching the job.
        processJobId = api.jobNew(processJobInputHash)
        log.debug "Launching job > ${processJobId}"

        // signal the new task status
        status = SUBMITTED
    }

    @Override
    void kill() {
        if( !processJobId ) { return }
        log.debug "Killing DnaNexus job with id: $processJobId"
        api.jobTerminate(processJobId)
    }

    @Override
    boolean checkIfRunning() {

        if( processJobId && isSubmitted() ) {

            def result = checkStatus()
            String state = result.state
            log.debug "Process ${task.name} > State: ${state}"

            if( state in ['idle', 'waiting_on_input', 'runnable', 'running', 'waiting_on_output','done','terminating'] ) {
                status = RUNNING
                return true
            }
        }

        return false
    }

    @Override
    boolean checkIfCompleted() {

        if( !isRunning() ) { return false }

        def result = checkStatus()
        String state = result.state
        if( !state ) { throw new IllegalStateException() }
        log.debug "Process ${task.name} > State: ${state}"

        if( !(state in ['done','terminating']) ) {
            return false
        }

        /*
         * Getting the exit code of the task's execution.
         */
        Integer exitCode = result.output?.exit_code
        if( exitCode != null ) {
            task.exitStatus = exitCode
            log.debug "Process ${task.name} > exit status > ${task.exitStatus}"
        }
        else {
            log.debug "Process ${task.name} > missing exit status"
        }

        /*
         * Getting the program output file.
         * When the 'echo' property is set, it prints out the task stdout
         */

        // the file that will receive the stdout
        Path taskOutputFile = task.getCmdOutputFile()
        if( !taskOutputFile.exists()) {
            log.warn "Process ${task.name} > output file does not exist: $taskOutputFile"
            task.stdout = '(none)'
        }
        else {
            log.debug "Process ${task.name} > out file: $taskOutputFile -- exists: ${taskOutputFile.exists()}; size: ${taskOutputFile.size()}"
            task.stdout = taskOutputFile
        }

        status = COMPLETED
        return true

    }


    private Map checkStatus() {

        long delta = System.currentTimeMillis() - lastStatusMillis
        if( delta < 15_000 ) {
            return lastStatusResult
        }

        if( !processJobId ) {
            throw new IllegalStateException()
        }

        def response = api.jobDescribe(processJobId)
        log.debug "Process ${task.name} > current result: ${response.toString()}\n"

        lastStatusMillis = System.currentTimeMillis()
        lastStatusResult = response
    }
}