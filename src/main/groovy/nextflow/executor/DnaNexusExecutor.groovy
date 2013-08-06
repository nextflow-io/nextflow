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
import com.dnanexus.DXAPI
import com.dnanexus.DXJSON
import com.fasterxml.jackson.databind.JsonNode
import com.fasterxml.jackson.databind.ObjectMapper
import com.fasterxml.jackson.databind.node.ObjectNode
import com.fasterxml.jackson.databind.node.ArrayNode
import groovy.util.logging.Slf4j
import nextflow.exception.MissingFileException
import nextflow.processor.TaskRun
import nextflow.util.DxFile
/**
 * Executes script.nf indicated in code.sh in the DnaNexus environment
 * -->  https://www.dnanexus.com/
 *
 * @author Beatriz Martin San Juan <bmsanjuan@gmail.com>
 */


@Slf4j
class DnaNexusExecutor extends AbstractExecutor {

    private static final COMMAND_OUT_FILENAME = '.command.out'

    private static final COMMAND_RUNNER_FILENAME = '.command.run'

    private static final COMMAND_ENV_FILENAME = '.command.env'

    private static final COMMAND_SCRIPT_FILENAME = '.command.sh'

    /*
     * Creation of Dnanexus link out of the id of an object.
     * @param objectId
     * @return ObjectNode
     */
    protected ObjectNode makeDXLink(String objectId) {
        return DXJSON.getObjectBuilder().put('$dnanexus_link', objectId).build();
    }


    /**
     * Creation of a link to a job and one of its fields out of the id of the first and
     * the required name of the second.
     * @param jobId
     * @param fieldName
     * @return ObjectNode
     */
    protected ObjectNode makeJbor(String jobId, String fieldName) {
        return DXJSON.getObjectBuilder().put("job", jobId).put("field", fieldName).build();
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
        def taskEnv = createEnvironmentString(task)
        log.debug "Creating Environment"


        /*
         * Saving the task script file in the appropriate task's working folder
         */
        File taskScript = new File(scratch, 'taskScript')
        taskScript.text = task.processor.normalizeScript(task.script.toString())
        log.debug "Creating script file for task > ${task.name}: ${taskScript.text}\n\n "


        /*
         * Uploading the task's script file.
         */
        Process scriptCmd = Runtime.getRuntime().exec("dx upload --brief ${taskScript.absolutePath}")
        BufferedReader uploadScript = new BufferedReader(new InputStreamReader(scriptCmd.getInputStream()))
        String scriptId = uploadScript.readLine().trim();
        log.debug "Uploading script file for task ${task.name} >> ${scriptId}"


        /*
         * Retrieving & uploading all the inputs already declared as parameters in the task.
         * Different method depending on the instance DxFile or File.
         */
        def map = task.code.delegate
        ObjectMapper mapper = new ObjectMapper();
        ArrayNode inputs = mapper.createArrayNode();

        map.each{ k, v ->
            if( v instanceof DxFile ) {
                inputs.add(makeDXLink(v.getId()));
                log.debug "Getting input DxFile ${k} for task ${task.name} >> Name: ${v} >> ${v.getId()}"
            }
            else if( v instanceof File ) {
                String path = String.valueOf(v)
                Process inputCmd = Runtime.getRuntime().exec("dx upload --brief ${path}")
                BufferedReader uploadInput = new BufferedReader(new InputStreamReader(inputCmd.getInputStream()))
                String inputId = uploadInput.readLine().trim();

                inputs.add(makeDXLink(inputId))
                log.debug "Uploading input file ${k} for task ${task.name} >> ${inputId}"
            }
            else {
                log.warn "Unsupported input type: $k --> $v"
            }
        }


        /*
         * Retrieving all the outputs already declared as parameters in the task.
         */
        ArrayNode outputs = mapper.createArrayNode()
        taskConfig.getOutputs().keySet().each { String name ->
             outputs.add(name)
        }


        /*
         * Building the ObjectNode which will be set in the job.
         * As parameters of this:
         *      - inputs --> String formed by all the names and ids of the inputs declared.
         *      - outputs --> String formed by all the names of the outputs declared.
         *      - taskname --> String with the name of the task.
         *      - taskScript --> Dnanexus link to the task's script file.
         */
        ObjectNode processJobInputHash = DXJSON.getObjectBuilder()
                .put("function", "process")
                .put("input", DXJSON.getObjectBuilder()
                        .put("inputs", inputs)
                        .put("outputs", outputs)
                        .put("taskName", task.name)
                        .put("taskScript", makeDXLink(scriptId))
                        .build())
                .build()
        //.put("taskEnv", '')
        log.debug "Creating job parameters"


        /*
         * Launching the job.
         */
        String processJobId = DXAPI.jobNew(processJobInputHash).get("id").textValue()
        log.debug "Launching job > ${processJobId}"
        task.jobId = processJobId


        /*
         * Waiting for the job to end while showing the state job's state and details.
         */
        JsonNode result = null
        String state = null
        log.debug "Waiting for the job"

        while( true ) {
            sleep( 15_000 )
            result = DXAPI.jobDescribe(processJobId)
            log.debug "Task ${task.name} -- current result: ${result.toString()}\n"

            state = result.get('state').textValue()
            if( state in ['idle', 'waiting_on_input', 'runnable', 'running', 'waiting_on_output', 'terminating'] ) {
                log.debug "State > ${state}"
                continue
            }
            log.debug "State > ${state}"
            break
        }


        /*
         * Getting the exit code of the task's execution.
         */
        String exitCode = result.get('output').get('exit_code').textValue()

        if( state == 'done' && exitCode?.isInteger() ) {
            task.exitCode = exitCode.toInteger()
            log.debug "Task's exit code > ${task.exitCode}"
        }

        /*
         * Getting the program output file.
         * When the 'echo' property is set, it prints out the task stdout
         */
        String outFileId = result.get('output')?.get('.command.out')?.textValue()
        if( !outFileId ) {
            log.warn "Unable to get task out file-id"
            task.output = '(unknown)'
        }
        else {
            log.debug "Downloading cmd out: $outFileId"
            def cmdOutFile = new File(scratch, '.command.out')
            DxHelper.downloadFile(outFileId, cmdOutFile)
            task.output = cmdOutFile

            if( taskConfig.echo ) {
                print cmdOutFile.text
            }
        }

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


    /**
     * Given the task and the name of one of the output files, it returns
     * all the files generated by the task's execution which matches the
     * file name.
     * The '-' stands for the script stdout, save to a file
     * @param task
     * @param fileName
     * @return  DxFile[] or DxFile
     */
    def collectResultFile( TaskRun task, String fileName ) {
        assert fileName
        assert task
        assert task.jobId

        if( fileName == '-' ) {
            return getStdOutFile(task)
        }

        JsonNode node = DXAPI.jobDescribe(task.jobId?.toString())
        def output = node.get('output')
        log.debug "Dx output: ${output.toString()}"

        return getFiles(output, fileName)
    }


    /**
     * Given the list JSON node containing the job output, return the {@code DxFile} instance
     * for the fine specified by the string {@code fileName}.
     * <p>
     *     When {@code fileName} contains one or more wildcards (star or question mark) a list of
     *     of files may be returned
     *
     * @param output
     * @param fileName
     * @return DxFile[] or DxFile
     */
    def getFiles( JsonNode outputs, String fileName ) {
        assert outputs != null
        assert fileName

        String filePattern = fileName.replace('?', '.?').replace('*', '.*')

        if( fileName == filePattern ) {
            String file = outputs.get(fileName)?.textValue()

            if( !file ) {
                throw new MissingFileException("Missing output file(s): '$fileName' expected by task: ${taskConfig.name}")
            }

            def result = new DxFile(id:file, name: fileName)
            log.debug "Result File >> ${result.getName()} : ${result.getId()}"

            return result
        }


        def result = []
        for( Map.Entry<String,JsonNode> entry : outputs.fields() ) {
            if( entry.key ==~/$filePattern/ ) {
                def fileId = entry.value?.textValue()
                log.debug "Result File >> ${fileName} >> ${entry.key} >> ${fileId}"
                def file = new DxFile(name: entry.key, id: entry.value?.textValue())
                result << file
            }
        }

        if( !result ) {
            throw new MissingFileException("Missing output file(s): '$fileName' expected by task: ${taskConfig.name}")
        }

        return result
    }
}