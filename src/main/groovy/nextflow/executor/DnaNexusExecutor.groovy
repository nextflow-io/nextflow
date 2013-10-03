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
import java.nio.file.Path

import com.dnanexus.DXAPI
import com.fasterxml.jackson.databind.JsonNode
import groovy.util.logging.Slf4j
import nextflow.fs.dx.DxPath
import nextflow.processor.TaskRun
import nextflow.util.DxHelper
/**
 * Executes script.nf indicated in code.sh in the DnaNexus environment
 * -->  https://www.dnanexus.com/
 *
 * @author Beatriz Martin San Juan <bmsanjuan@gmail.com>
 */


@Slf4j
class DnaNexusExecutor extends AbstractExecutor {

    private static final COMMAND_OUT_FILENAME = '.command.out'

    private static final COMMAND_IN_FILENAME = '.command.in'

    private static final COMMAND_ENV_FILENAME = '.command.env'

    private static final COMMAND_SCRIPT_FILENAME = '.command.sh'


    private List getInputFiles(TaskRun task) {

        def map = task.code.delegate
        def inputs = []

        map.each{ k, v ->
            if( v instanceof DxPath ) {
                log.debug "Getting input DxPath ${k} for task ${task.name} >> ${v} >> ${v.getFileId()}"
                inputs.add( (v as DxPath).getFileId() );
            }
        }

        return inputs
    }

    private List getOutputFiles( TaskRun task ) {

        new ArrayList(taskConfig.getOutputs().keySet())

    }


    /**
     * Launches the task
     * @param task
     */
    @Override
    void launchTask(TaskRun task) {

        /*
         * Setting the work directory
         */
        final scratch = task.workDirectory
        log.debug "Lauching task > ${task.name} -- work folder: $scratch"


        /*
         * Saving the environment to a file.
         */
        def taskEnvFile = scratch.resolve(COMMAND_ENV_FILENAME);
        createEnvironmentFile(task, taskEnvFile)

        /*
         * In case there's a task input file.
         */
        Path taskInputFile = null
        if (task.input) {

            /*
             * Saving the task input file in the appropriate task's working folder
             */
            taskInputFile = scratch.resolve(COMMAND_IN_FILENAME)
            taskInputFile.text = task.input

        }


        /*
         * Saving the task script file in the appropriate task's working folder
         */
        Path taskScriptFile = scratch.resolve(COMMAND_SCRIPT_FILENAME)
        taskScriptFile.text = task.processor.normalizeScript(task.script.toString())
        log.debug "Creating script file for task > ${task.name}\n\n "


        // input job params
        def obj = [:]
        obj.task_name = task.name
        obj.task_script = (taskScriptFile as DxPath).getFileId()
        //obj.task_env = (taskEnvFile as DxPath).getFileId()
        if( taskInputFile ) {
            obj.task_input = (taskInputFile as DxPath).getFileId()
        }
        obj.input_files = getInputFiles(task)
        obj.output_files = getOutputFiles(task)

        // create the input parameters for the job to be executed
        def processJobInputHash = createInputObject( obj, (String)taskConfig.instaceType )
        log.debug "New job parameters: ${processJobInputHash}"

        /*
         * Launching the job.
         */
        String processJobId = DXAPI.jobNew(processJobInputHash).get("id").textValue()
        log.debug "Launching job > ${processJobId}"
        task.jobId = processJobId


        /*
         * Waiting for the job to end while showing the state job's state and details.
         */
        log.debug "Waiting for the job > $processJobId"
        Map result = waitForJobResult(task)

        /*
         * Getting the exit code of the task's execution.
         */
        Integer exitCode = result.output?.exit_code
        if( result.state == 'done' && exitCode != null ) {
            task.exitCode = exitCode
            log.debug "Task exit code > ${task.exitCode}"
        }


        /*
         * Getting the program output file.
         * When the 'echo' property is set, it prints out the task stdout
         */
        // the file that will receive the stdout
        Path taskOutputFile = scratch.resolve(COMMAND_OUT_FILENAME)
        if( !taskOutputFile.exists()) {
            log.warn "Task output file does not exist: $taskOutputFile"
            task.output = '(unknown)'
        }
        else {
            log.debug "Task out file: $taskOutputFile -- exists: ${taskOutputFile.exists()}; size: ${taskOutputFile.size()}\n ${taskOutputFile.text} "
            task.output = taskOutputFile

            if( taskConfig.echo ) {
                print taskOutputFile.text
            }
        }
    }

    def Map waitForJobResult( TaskRun task ) {
        JsonNode result;

        while( true ) {
            sleep( 15_000 )
            result = DXAPI.jobDescribe(task.jobId as String)
            log.debug "Task ${task.name} -- current result: ${result.toString()}\n"

            String state = result.get('state').textValue()
            if( state in ['idle', 'waiting_on_input', 'runnable', 'running', 'waiting_on_output', 'terminating'] ) {
                log.debug "State > ${state}"
                continue
            }
            log.debug "State > ${state}"
            break
        }

        return DxHelper.jsonToObj(result)
    }


    def resolveInputFile( def obj ) {
        if( obj instanceof DxPath ) {
            // call getFileId() to force resolving it to te real file before get as simple path name
            obj.getFileId()
            // now return the the file name as path
            return obj.getFileName()
        }
        else {
            obj
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

    /**
     * Returns the output of the task.
     * @param task
     * @return task.output
     */
    @Override
    def getStdOutFile(TaskRun task) {
        task.@output
    }


}