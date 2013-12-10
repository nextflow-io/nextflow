package nextflow.processor
import java.nio.file.Path
import java.util.concurrent.atomic.AtomicInteger

import groovy.transform.InheritConstructors
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.operator.DataflowEventAdapter
import groovyx.gpars.dataflow.operator.DataflowOperator
import groovyx.gpars.dataflow.operator.DataflowProcessor
import groovyx.gpars.dataflow.operator.PoisonPill
import nextflow.util.CacheHelper
import nextflow.util.FileHelper
/**
 * Defines the 'merge' processing policy
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@Slf4j
@InheritConstructors
class MergeTaskProcessor extends TaskProcessor {

    /*
     * Collect the list of hashCode of each {@code TaskRun} object
     * required by this merge operation
     */
    private List<Integer> mergeHashesList

    /*
     * The folder where the merge will be execute
     */
    private Path mergeTempFolder

    /*
     * Count the total number of {@code TaskRun} operation
     * executed by the final merge operation
     */
    private AtomicInteger mergeIndex = new AtomicInteger()

    /*
     * The script which will execute the merge operation.
     * This script is created by collecting all the single {@code TaskRun} scripts serially
     */
    private def mergeScript = new StringBuilder()

    /* Collect all the input values used by executing this task */
    private Map<InParam,List> inputsCollector = [:]

    /*
     * The merge task is composed by two operator, the first creates a single 'script' to be executed by second one
     */
    @Override
    protected void createOperator() {
        log.debug "Starting merge > ${name}"

        def wrapper = createCallbackWrapper(taskConfig.inputs.size(), this.&mergeScriptCollector)
        mergeHashesList = new LinkedList<>()
        mergeTempFolder = FileHelper.createTempFolder(session.workDir)

        // initialize the output values collector
        taskConfig.inputs.each { InParam param ->
            inputsCollector[param] = []
        }


        def params = [inputs: new ArrayList(taskConfig.inputs.channels), outputs: new ArrayList(taskConfig.outputs.channels), listeners: [new MergeProcessorInterceptor()] ]
        processor = new DataflowOperator(group, params, wrapper)
        session.allProcessors.add(processor)

        // -- start it
        processor.start()
    }

    protected void mergeTaskRun(TaskRun task) {

        task.stagedProvider = this.&stagedProvider

        /*
         * initialize the *only* outputs for this task instance
         * (inputs have been managed during the scripts collection stage)
         */
        taskConfig.outputs.each { OutParam param -> task.setOutput(param) }

        // -- create the unique hash number for this tasks,
        //    collecting the id of all the executed runs
        //    and sorting them
        def hasher = CacheHelper.hasher(session.uniqueId)
        mergeHashesList.sort()
        mergeHashesList.each { Integer entry ->  hasher = CacheHelper.hasher(hasher,entry) }
        def hash = hasher.hash()
        log.trace "Merging task > $name -- hash: $hash"

        Path folder = FileHelper.getWorkFolder(session.workDir, hash)
        log.trace "Merging task > $name -- trying cached: $folder"

        def cached = session.cacheable && taskConfig.cacheable && checkCachedOutput(task,folder)
        if( !cached ) {

            folder = createTaskFolder(folder, hash)
            log.info "Running merge > ${name}"

            // -- set the folder where execute the script
            task.workDirectory = folder

            // -- set the aggregate script to be executed
            task.script = mergeScript.toString()

            // -- submit task for execution !
            submitTask( task )

        }

    }

    protected void mergeScriptCollector( List values ) {
        final currentIndex = mergeIndex.incrementAndGet()
        log.info "Collecting task > ${name} ($currentIndex)"

        // -- map the inputs to a map and use to delegate closure values interpolation
        def keys = []
        def stdin = null
        def contextMap = new DelegateMap(ownerScript)
        Map<FileInParam,List<FileHolder>> filesMap = [:]
        Map<String,String> environment = [:]
        int count = 0

        taskConfig.inputs.eachWithIndex { InParam param, int index ->

            final val = values[index]

            // define the *context* against which the script will be evaluated
            switch( param ) {
                case ValueInParam:
                    contextMap[param.name] = val
                    break

                case StdInParam:
                    stdin = val
                    break

                case FileInParam:
                    // all the files to be staged
                    def fileParam = (FileInParam)param
                    def normalized = normalizeInputToFiles(val,count)
                    def resolved = expandWildcards( fileParam.filePattern, normalized )
                    filesMap[fileParam] = resolved
                    count += resolved.size()
                    // set the context
                    contextMap[param.name] = singleItemOrList( resolved )
                    // set to *val* so that the list is added to the map of all inputs
                    val = resolved
                    break

                case EnvInParam:
                    // the environment variables for this 'iteration'
                    environment[param.name] = val
                    break

                default:
                    log.debug "Task $name > unknown input param type: ${param?.class?.simpleName}"
            }

            // store all the inputs
            inputsCollector.get(param).add( val )

            // add all the input name-value pairs to the key generator
            keys << param.name << val
        }


        /*
         * initialize the task code to be executed
         */
        Closure scriptClosure = this.code.clone() as Closure
        scriptClosure.delegate = contextMap
        scriptClosure.setResolveStrategy(Closure.DELEGATE_FIRST)

        def commandToRun = normalizeScript(scriptClosure.call()?.toString())
        def interpreter = fetchInterpreter(commandToRun)


        /*
         * create a unique hash-code for this task run and save it into a list
         * which maintains all the hashes for executions making-up this merge task
         */
        keys << commandToRun << 7
        mergeHashesList << CacheHelper.hasher(keys).hash().asInt()

        // section marker
        mergeScript << "# task '$name' ($currentIndex)" << '\n'

        // add the files to staged
        if( filesMap ) {
            mergeScript << executor.stagingFilesScript(filesMap)
        }

        // add the variables to be exported
        if( environment ) {
            mergeScript << bashEnvironmentScript(environment)
        }

        /*
         * save the script to execute into a separate unique-named file
         */
        final index = currentIndex
        final scriptName = ".merge_command.sh.${index.toString().padLeft(4,'0')}"
        final scriptFile = mergeTempFolder.resolve(scriptName)
        scriptFile.text = commandToRun

        // the command to launch this command
        def scriptCommand = scriptName

        // check if some input have to be send
        if( stdin ) {
            final inputName = ".merge_command.input.$index"
            final inputFile = mergeTempFolder.resolve( inputName )
            inputFile.text = stdin

            // pipe the user input to the user command
            scriptCommand = "$scriptCommand < ${inputName}"

            // stage the input file
            mergeScript << executor.stageInputFileScript(inputFile, inputName) << '\n'
        }

        // stage this script itself
        mergeScript << executor.stageInputFileScript(scriptFile, scriptName) << '\n'

        // create a unique script collecting all the commands
        mergeScript << interpreter << ' ' << scriptCommand << '\n'

    }

    protected Map<FileInParam,List<FileHolder>> stagedProvider() {
        (Map<FileInParam,List<FileHolder>>) inputsCollector.findAll { it.key instanceof FileInParam }
    }



    /**
     * A task of type 'merge' binds the output when it terminates it's work, i.e. when
     * it receives a 'poison pill message that will stop it
     */
    class MergeProcessorInterceptor extends DataflowEventAdapter {


        @Override
        public List<Object> beforeRun(final DataflowProcessor processor, final List<Object> messages) {
            log.trace "Before run > ${name} -- messages: ${messages}"
            return messages;
        }

        @Override
        void afterRun(DataflowProcessor processor, List<Object> MOCK_MESSAGES) {
            log.trace "After run > ${name}"
        }

        @Override
        public Object messageArrived(final DataflowProcessor processor, final DataflowReadChannel<Object> channel, final int index, final Object message) {
            log.trace "Received message > task: ${name}; channel: $index; value: $message"
            return message
        }

        public Object controlMessageArrived(final DataflowProcessor processor, final DataflowReadChannel<Object> channel, final int index, final Object message) {
            if( log.isTraceEnabled() ) {
                def inputs = taskConfig.inputs.names
                log.trace "Control message arrived > ${name} -- ${inputs[index]} => ${message}"
            }

            // intercepts 'poison pill' message e.g. termination control message
            // in order the launch the task execution on the underlying system
            if( message == PoisonPill.instance )  {

                if( mergeIndex.get()>0 ) {
                    receivedPoisonPill = true
                    instanceCount.incrementAndGet()
                    def task = MergeTaskProcessor.this.createTaskRun()
                    try {
                        MergeTaskProcessor.this.mergeTaskRun(task)
                    }
                    catch( Throwable e ) {
                        handleException(e, task)
                    }

                    return StopQuietly.instance
                }

                log.warn "No data collected by task > $name -- Won't execute it. Something may be wrong in your execution flow"

            }

            return message
        }

        public boolean onException(final DataflowProcessor processor, final Throwable e) {
            handleException(e)
        }

    }


}
