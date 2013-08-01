package nextflow.executor
import com.dnanexus.DXAPI
import com.dnanexus.DXJSON
import com.fasterxml.jackson.databind.JsonNode
import com.fasterxml.jackson.databind.node.ObjectNode
import groovy.util.logging.Slf4j
import nextflow.exception.MissingFileException
import nextflow.processor.TaskRun

/**
 * Created with IntelliJ IDEA.
 * User: bmartin
 * Date: 7/18/13
 * Time: 1:54 PM
 * To change this template use File | Settings | File Templates.
 */

@Slf4j
class DnaNexusExecutor extends AbstractExecutor {

    private static final COMMAND_OUT_FILENAME = '.command.out'

    private static final COMMAND_RUNNER_FILENAME = '.command.run'

    private static final COMMAND_ENV_FILENAME = '.command.env'

    private static final COMMAND_SCRIPT_FILENAME = '.command.sh'


    protected ObjectNode makeDXLink(String objectId) {
        return DXJSON.getObjectBuilder().put("\$dnanexus_link", objectId).build();
    }

    protected ObjectNode makeJbor(String jobId, String fieldName) {
        return DXJSON.getObjectBuilder().put("job", jobId).put("field", fieldName).build();
    }


    @Override
    void launchTask(TaskRun task) {
        //To change body of implemented methods use File | Settings | File Templates.

        final scratch = task.workDirectory
        log.debug "Lauching task > ${task.name} -- work folder: $scratch"


        /*
         * save the environment to a file
         */
        def taskEnv = createEnvironmentString(task)
        log.debug "Creating Environment"


        /*
         * Saving the main script file
         */
        File taskScript = new File('taskScript')
        taskScript.text = task.processor.normalizeScript(task.script.toString())
        println("PATH: ${taskScript.absolutePath}")
        log.debug "Creating script file > ${taskScript.name}"


        /*
         * Uploading the main script file
         */
        Process scriptCmd = Runtime.getRuntime().exec("dx upload --brief ${taskScript.absolutePath}")
        BufferedReader uploadScript = new BufferedReader(new InputStreamReader(scriptCmd.getInputStream()))
        String scriptId = uploadScript.readLine().trim(); // In general, there'd be one line per uploaded file
        log.debug "Uploading script file >> ${scriptId}"


        /**
         * Retrieving & uploading all the inputs of the task
         */
        def map = task.code.delegate
        String inputs = "( "

        map.each{ k, v ->
            String path = String.valueOf(v)

            if( v instanceof DxFile ) {
                String link = makeDXLink(v.getId())

                inputs = inputs + "[${k}]=\'${link}\' "
                log.debug "Getting input DxFile ${k} >> ${v.getId()}"
            }
            else if( v instanceof File ) {
                Process inputCmd = Runtime.getRuntime().exec("dx upload --brief ${path}")
                BufferedReader uploadInput = new BufferedReader(new InputStreamReader(inputCmd.getInputStream()))
                String inputId = uploadInput.readLine().trim();
                String link = makeDXLink(inputId)

                inputs = inputs + "[${k}]=\'${link}\' "
                log.debug "Uploading input file ${k} >> ${inputId}"
            }
            else {
                log.warn "Unsupported input type: $k --> $v"
            }
        }

        inputs = inputs + ")"

        String outputs = "( "

        taskConfig.getOutputs().keySet().each { String name ->
             outputs= outputs + "\"${name}\" "
        }
        outputs = outputs + ")"


        /*
         * Creating the job with a input's Map
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
         * Launching the job
         */
        String processJobId = DXAPI.jobNew(processJobInputHash).get("id").textValue()
        log.debug "Launching job > ${processJobId}"


        /*
         * Waiting for the job to end
         */
        JsonNode result = null
        String state = null
        log.debug "Waiting for the job"

        while( true ) {
            sleep( 10_000 )
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
         * Getting the task's outputs
         */
        String exitCode = result.get('output').get('exit_code').textValue()

        if( state == 'done' && exitCode?.isInteger() ) {
            task.exitCode = exitCode.toInteger()
            log.debug "Task's exit code > ${task.exitCode}"
        }

        task.output=null

        /*
         * Fake job exit status -- this should come from the job executed in the cloud
         */
        task.exitCode = 0


        /*
         * Fake job stdout -- this have to be download from dna-nexus
         */

        File fileOut = new File(scratch, COMMAND_OUT_FILENAME)
        fileOut.text = '(fake result)\n\n'
        task.output = fileOut

        /*
         * Show the job result when the 'echo' flag is set
         * This have to be replaced by the API call /job-xxx/streamLog -- http://goo.gl/3EFcZW
         */
        if( taskConfig.echo ) {
            System.out.print( fileOut.text )
        }


    }


    @Override
    def getStdOutFile(TaskRun task) {
        task.@output
    }


    def collectResultFile( TaskRun task, String fileName ) {
        assert fileName
        assert task
        assert task.jobId

        // the '-' stands for the script stdout, save to a file
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
     */
    def getFiles( JsonNode outputs, String fileName ) {

        if(fileName.contains('*') || fileName.contains('?')){
            String output_string = outputs.toString()
            String[] array = output_string.substring(1,output_string.length()-2).split(",")

            for(int i=0; i< array.length; i++){
                array[i] = array[i].substring(0,array[i].indexOf(":")).replace("\"","")
            }

            def result = []
            String expression

            if(fileName.contains('*'))
                expression =  fileName.replace(".","\\.").replace("*",".*")
            else if (fileName.contains('?'))
                expression =  fileName.replace(".","\\.").replace("?",".?")


            for(int i=0; i<array.length; i++ ) {
                if(array[i] =~ /^$expression$/){
                    String fileId = outputs.get(array[i])?.textValue()
                    DxFile file = new DxFile(id:fileId, name: array[i])
                    result.add(file)
                    log.debug "Result File >> ${fileName} >> ${array[i]} >> ${fileId}"
                }
            }

            if( !result ) {
                throw new MissingFileException("Missing output file(s): '$fileName' expected by task: ${taskConfig.name}")
            }

            return result
        }
        else{
            String file = outputs.get(fileName)?.textValue()

            if( !file ) {
                throw new MissingFileException("Missing output file(s): '$fileName' expected by task: ${taskConfig.name}")
            }

            def result = new DxFile(id:file, name: fileName)
            log.debug "Result File >> ${result.getName()} : ${result.getId()}"

            return result
        }
    }
}