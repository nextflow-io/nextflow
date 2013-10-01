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
import com.dnanexus.DXJSON
import com.fasterxml.jackson.databind.JsonNode
import com.fasterxml.jackson.databind.node.ObjectNode
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

    /*
     * Creation of Dnanexus link out of the id of an object.
     * @param objectId
     * @return ObjectNode
     */
    protected static ObjectNode makeDXLink(String objectId) {
        return DXJSON.getObjectBuilder().put('$dnanexus_link', objectId).build();
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
        obj.output_files = new ArrayList(taskConfig.getOutputs().keySet())

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

//        if(!taskInputId.equals(null)){
//            root.input.taskInput = makeDXLink(taskInputId)
//        }
//        root.input.taskScript = makeDXLink(scriptId)
//        root.input.taskEnv = env
//        root.input.taskName = taskName
//        root.input.outputs = outputs
//        root.input.inputs = inputObj

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


//    /**
//     * Given the task and the name of one of the output files, it returns
//     * all the files generated by the task's execution which matches the
//     * file name.
//     * The '-' stands for the script stdout, save to a file
//     * @param task
//     * @param fileName
//     * @return  DxFile[] or DxFile
//     */
//    def collectResultFile( TaskRun task, String fileName ) {
//        assert fileName
//        assert task
//        assert task.jobId
//
//        if( fileName == '-' ) {
//            return getStdOutFile(task)
//        }
//
//        JsonNode node = DXAPI.jobDescribe(task.jobId?.toString())
//        def output = node.get('output')
//        log.debug "Dx output: ${output.toString()}"
//
//        return getFiles(output, fileName)
//    }
//
//
//    /**
//     * Given the list JSON node containing the job output, return the {@code DxFile} instance
//     * for the fine specified by the string {@code fileName}.
//     * <p>
//     *     When {@code fileName} contains one or more wildcards (star or question mark) a list of
//     *     of files may be returned
//     *
//     * @param output
//     * @param fileName
//     * @return DxFile[] or DxFile
//     */
//    def getFiles( JsonNode outputs, String fileName ) {
//        assert outputs != null
//        assert fileName
//
//        String filePattern = fileName.replace('?', '.?').replace('*', '.*')
//
//        if( fileName == filePattern ) {
//            String file = outputs.get(fileName)?.textValue()
//
//            if( !file ) {
//                throw new MissingFileException("Missing output file(s): '$fileName' expected by task: ${taskConfig.name}")
//            }
//
//            def result = new DxFile(id:file, name: fileName)
//            log.debug "Result File >> ${result.getName()} : ${result.getId()}"
//
//            return result
//        }
//
//
//        def result = []
//        for( Map.Entry<String,JsonNode> entry : outputs.fields() ) {
//            if( entry.key ==~/$filePattern/ ) {
//                def fileId = entry.value?.textValue()
//                log.debug "Result File >> ${fileName} >> ${entry.key} >> ${fileId}"
//                def file = new DxFile(name: entry.key, id: entry.value?.textValue())
//                result << file
//            }
//        }
//
//        if( !result ) {
//            throw new MissingFileException("Missing output file(s): '$fileName' expected by task: ${taskConfig.name}")
//        }
//
//        return result
//    }
}