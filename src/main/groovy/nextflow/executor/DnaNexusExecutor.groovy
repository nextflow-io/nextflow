package nextflow.executor
import com.dnanexus.DXJSON
import com.dnanexus.DXAPI
import com.fasterxml.jackson.databind.node.ObjectNode
import com.fasterxml.jackson.databind.JsonNode
import com.fasterxml.jackson.databind.ObjectMapper
import groovy.io.FileType
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
        println("OUTPUTS: ${outputs}")

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

    /*    ObjectNode jobOutput = DXJSON.getObjectBuilder()
                .put("output_file", makeJbor(processJobId, "file1"))
                .build();

        ObjectMapper mapper = new ObjectMapper();
        File final_output= new File("job_output.json")
        mapper.writeValue(final_output, jobOutput)
        log.debug "Final output >> ${final_output.text}"


        task.output = final_output                  */
        task.jobId = processJobId

        /*
        * Show the job result when the 'echo' flag is set
        * This have to be replaced by the API call /job-xxx/streamLog -- http://goo.gl/3EFcZW
        */
    /*    if( taskConfig.echo ) {
            println( final_output.toString()) //fileOut.text )
        }*/
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

        JsonNode node = DXAPI.jobDescribe(task.jobId)

        if(fileName.startsWith('*.')){
            String outputs = node.get('output')
            String[] array = outputs.substring(1,outputs.length()-2).split(",")

            for(int i=0; i< array.length; i++){
                array[i] = array[i].substring(0,array[i].indexOf(":")).replace("\"","")
            }

            String name
            String expression =  fileName.replace(".","\\.").replace("*",".*")

            for(int i=0; i<array.length; i++ ) {
                if(array[i] =~ /$expression/){
                    name = array[i]
                }
            }

            String fileFinal = node.get('output').get(name).textValue()
            def result = new DxFile(id:fileFinal, name: name)
            log.debug "Result File >> ${fileName} >> ${name} >> ${fileFinal}"

            return result

        }else{
            String file = node.get('output').get(fileName).textValue()

            def result = new DxFile(id:file, name: fileName)
            log.debug "Result File >> ${result.getName()} : ${result.getId()}"

            return result
        }
    }
}