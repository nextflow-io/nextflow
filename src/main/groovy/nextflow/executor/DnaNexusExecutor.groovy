/*
* Copyright (c) 2012, the authors.
*
* This file is part of 'Nextflow'.
*
* Nextflow is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* Nextflow is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with Nextflow. If not, see <http://www.gnu.org/licenses/>.
*/

package nextflow.executor
import java.nio.file.Files
import java.nio.file.Path

import com.dnanexus.DXAPI
import com.fasterxml.jackson.databind.JsonNode
import groovy.util.logging.Slf4j
import nextflow.fs.dx.DxPath
import nextflow.processor.FileInParam
import nextflow.processor.FileOutParam
import nextflow.processor.TaskConfig
import nextflow.processor.TaskHandler
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import nextflow.util.DxHelper

/**
 * Executes script.nf indicated in dxapp.sh in the DnaNexus environment
 *
 * See https://www.dnanexus.com/
 *
 * @author Beatriz Martin San Juan <bmsanjuan@gmail.com>
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */


@Slf4j
class DnaNexusExecutor extends AbstractExecutor {

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
        log.debug "Lauching task > ${task.name} -- work folder: $scratch"

        /*
         * Saving the environment to a file.
         */
        Map environment = task.processor.getProcessEnvironment()
        environment.putAll( task.getInputEnvironment() )
        final taskEnvFile = task.getCmdEnvironmentFile()
        taskEnvFile.text = TaskProcessor.bashEnvironmentScript(environment)

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
        log.debug "Creating script file for task > ${task.name}\n\n "

        /*
         * create the stage inputs script
         */
        def inputFiles = task.getInputsByType(FileInParam)
        def stageScript = stagingFilesScript( inputFiles, '; ' )

        /*
         * the DnaNexus job params object
         */
        def obj = [:]
        obj.task_name = task.name
        obj.task_script = (taskScriptFile as DxPath).getFileId()
        obj.task_env = (taskEnvFile as DxPath).getFileId()
        if( taskInputFile ) {
            obj.task_input = (taskInputFile as DxPath).getFileId()
        }
        obj.stage_inputs = stageScript
        obj.output_files = taskConfig.getOutputs().ofType(FileOutParam).collect { it.getName() }

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

    def static JsonNode createInputObject( Map inputObj, String instanceType ) {

        def root = [:]

        if(instanceType){
            def process = [ instanceType: instanceType ]
            root.systemRequirements = [process: process]
        }

        root.input = inputObj
        root.function = "process"

        return DxHelper.objToJson(root)

    }

}

/**
 * Handle a job execution in the DnaNexus platform
 */
@Slf4j
class DxTaskHandler extends TaskHandler {

    final Map inputParams

    final DnaNexusExecutor executor

    private String processJobId

    private long lastStatusMillis

    private Map lastStatusResult


    protected DxTaskHandler(TaskRun task, TaskConfig config, DnaNexusExecutor executor, Map params) {
        super(task, config)
        this.taskConfig = config
        this.inputParams = params
        this.executor = executor
    }


    @Override
    void submit() {

        // create the input parameters for the job to be executed
        def processJobInputHash = executor.createInputObject( inputParams, (String)taskConfig.instaceType )
        log.debug "New job parameters: ${processJobInputHash}"

        // Launching the job.
        processJobId = DXAPI.jobNew(processJobInputHash).get("id").textValue()
        log.debug "Launching job > ${processJobId}"

    }

    @Override
    void kill() {
        if( !processJobId ) { return }
        log.debug "Killing DnaNexus job with id: $processJobId"
        DXAPI.jobTerminate(processJobId)
    }

    @Override
    boolean checkIfRunning() {

        if( !isNew() ) {
            return true
        }

        def result = checkStatus()
        String state = result.state
        log.debug "Task ${task.name} > State: ${state}"

        if( state in ['idle', 'waiting_on_input', 'runnable', 'running', 'waiting_on_output'] ) {
            status = Status.RUNNING
            return true
        }

        return false
    }

    @Override
    boolean checkIfTerminated() {

        if( isTerminated() ) { return true }

        if( !isRunning() ) { return false }

        def result = checkStatus()
        String state = result.state
        if( !state ) { throw new IllegalStateException() }
        log.debug "Task ${task.name} > State: ${state}"

        if( !(state in ['done','terminating']) ) {
            return false
        }

        /*
         * Getting the exit code of the task's execution.
         */
        Integer exitCode = result.output?.exit_code
        if( exitCode != null ) {
            task.exitCode = exitCode
            log.debug "Task ${task.name} > exit code > ${task.exitCode}"
        }
        else {
            log.debug "Task ${task.name} > missing exit code"
        }

        /*
         * Getting the program output file.
         * When the 'echo' property is set, it prints out the task stdout
         */

        // the file that will receive the stdout
        Path taskOutputFile = task.getCmdOutputFile()
        if( !taskOutputFile.exists()) {
            log.warn "Task ${task.name} > output file does not exist: $taskOutputFile"
            task.stdout = '(none)'
        }
        else {
            log.debug "Task ${task.name} > out file: $taskOutputFile -- exists: ${taskOutputFile.exists()}; size: ${taskOutputFile.size()}\n ${taskOutputFile.text} "
            task.stdout = taskOutputFile
        }

        status = Status.TERMINATED
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

        def response = DXAPI.jobDescribe(processJobId)
        log.debug "Task ${task.name} > current result: ${response.toString()}\n"

        lastStatusMillis = System.currentTimeMillis()
        lastStatusResult = DxHelper.jsonToObj(response)
    }
}