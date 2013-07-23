package nextflow.executor
import com.dnanexus.DXJSON
import com.fasterxml.jackson.databind.node.ObjectNode
import groovy.util.logging.Slf4j
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

        /*
         * save the main script file
         */
        File taskScript = new File('taskScript')
        taskScript.text = task.processor.normalizeScript(task.script.toString())
        String path= taskScript.absolutePath

        Process uploadCmd = Runtime.getRuntime().exec("dx upload --brief ${path}")
        // Can specify multiple files above

        BufferedReader uploadOutput = new BufferedReader(new InputStreamReader(uploadCmd.getInputStream()))
        String fileId = uploadOutput.readLine().trim(); // In general, there'd be one line per uploaded file

        log.debug "File >> ${fileId}"
        println("END OF PROGRAM")


       // ObjectMapper mapper = new ObjectMapper();
    /*    ObjectNode processJobInputHash = DXJSON.getObjectBuilder()
                .put("function", "process")
                .put("input", DXJSON.getObjectBuilder()
                                    .put("taskName", task.name)
                                    .put("taskScript", makeDXLink(fileId))
                                    .build())
                .build()            **/
        //.put("taskEnv", '')


        /*
         * launch the job
         */
       // String processJobId = DXAPI.jobNew(processJobInputHash).get("id").textValue()


        /*
         * Wait for task execution
         */


     /*   JsonNode result = null
        String state = null
        while( true ) {
            sleep( 10_000 )
            result = DXAPI.jobDescribe(processJobId)
            log.debug "Task ${task.name} -- current result: ${result.toString()}\n"

            state = result.get('state').textValue()
            if( state in ['idle', 'waiting_on_input', 'runnable', 'running', 'waiting_on_output', 'terminating'] ) {
                continue
            }
            break
        }

        String exitCode = result.get('output').get('exit_code').textValue()
        log.debug "task exit code as string: $exitCode"

        if( state == 'done' && exitCode?.isInteger() ) {
            task.exitCode = exitCode.toInteger()
        }       */


//        def output = result.get('output').textValue()
//        def failureReason = result.get('failureReason').textValue()
//        def failureMessage = result.get('failureMessage').textValue()
//
//        log.debug "task terminated > output: $output -- failureReason: $failureReason -- failureMessage: $failureMessage"


     //   DXAPI.fileDownload('xx')


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
}
