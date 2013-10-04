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
import nextflow.exception.InvalidExitException
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

    private List<Integer> mergeHashesList

    private Path mergeTempFolder

    private AtomicInteger mergeIndex = new AtomicInteger()

    private def mergeScript = new StringBuilder()

    private Map<String,List> valuesCollector = [:]

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
        taskConfig.outputs.ofType(ValueOutParam).each { ValueOutParam param ->
            if( taskConfig.inputs.ofType(ValueInParam).find { it.name == param.name } ) {
                valuesCollector[param.name] = []
            }
            else {
                log.warn "Not a valid output parameter: '${param.name}' -- only values declared as input can be used in the output section"
            }
        }

        def params = [inputs: new ArrayList(taskConfig.inputs.channels), outputs: new ArrayList(taskConfig.outputs.channels), listeners: [new MergeProcessorInterceptor()] ]
        processor = new DataflowOperator(group, params, wrapper)
        session.allProcessors.add(processor)

        // -- start it
        processor.start()
    }


    protected void mergeScriptCollector( List values ) {
        final currentIndex = mergeIndex.incrementAndGet()
        log.info "Collecting task > ${name} ($currentIndex)"

        // -- map the inputs to a map and use to delegate closure values interpolation
        def keys = []
        def stdin = null
        def contextMap = new DelegateMap(ownerScript)
        Map<FileInParam,Object> filesMap = [:]
        Map<String,String> environment = [:]

        taskConfig.inputs.eachWithIndex { InParam param, int index ->

            // define the *context* against which the script will be evaluated
            if( param instanceof ValueInParam ) {
                contextMap[param.name] = values[index]
            }
            // define the *stdin* text
            else if( param instanceof StdInParam ) {
                stdin = values[index]
            }
            // all the files to be staged
            else if( param instanceof FileInParam ) {
                filesMap[param] = values[index]
            }
            // the environment variables for this 'iteration'
            else if( param instanceof EnvInParam ) {
                environment[param.name] = values[index]
            }
            else {
                log.debug "Task $name > unknown input param type: ${param?.class?.simpleName}"
            }

            // add all the input name-value pairs to the key generator
            keys << param.name << values[index]
        }

        /*
         * initialize the task code to be executed
         */
        Closure scriptClosure = this.code.clone() as Closure
        scriptClosure.delegate = contextMap
        scriptClosure.setResolveStrategy(Closure.DELEGATE_FIRST)

        def commandToRun = normalizeScript(scriptClosure.call()?.toString())
        def interpreter = fetchInterpreter(commandToRun)

        // collect the output values
        valuesCollector.keySet().each { key ->
            if( contextMap.containsKey(key) ) {
                valuesCollector[key] << contextMap[key]
            }
        }


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
            mergeScript << executor.stagingFilesScript( filesMap )
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

    protected void mergeTaskRun(TaskRun task) {

        try {

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

            Path folder = FileHelper.createWorkFolder(session.workDir, hash)
            log.trace "Merging task > $name -- trying cached: $folder"

            def cached = session.cacheable && taskConfig.cacheable && checkCachedOutput(task,folder)
            if( !cached ) {

                folder = createTaskFolder(folder, hash)
                log.info "Running merge > ${name}"

                // -- set the folder where execute the script
                task.workDirectory = folder

                // -- set the aggregate script to be executed
                task.script = mergeScript.toString()

                // -- run it !
                launchTask( task )

                // -- save the exit code
                folder.resolve('.exitcode').text = task.exitCode

                // -- check if terminated successfully
                boolean success = (task.exitCode in taskConfig.validExitCodes)
                if ( !success ) {
                    throw new InvalidExitException("Task '${task.name}' terminated with an invalid exit code: ${task.exitCode}")
                }

                collectOutputs(task)
                bindOutputs(task)
            }
        }
        finally {
            task.status = TaskRun.Status.TERMINATED
        }
    }

    @Override
    protected void collectOutputs( TaskRun task, OutParam param ) {

        if( param instanceof ValueOutParam ) {
            task.setOutput(param, valuesCollector[param.name])
            return
        }

        // fallback on the default behavior
        super.collectOutputs(task, param)
    }




    /**
     * A task of type 'merge' binds the output when it terminates it's work, i.e. when
     * it receives a 'poison pill message that will stop it
     */
    class MergeProcessorInterceptor extends DataflowEventAdapter {

        @Override
        public void afterStart(final DataflowProcessor processor) {
            // increment the session sync
            def val = session.taskRegister()
            log.trace "After start > register phaser '$name' :: $val"
        }

        @Override
        void afterRun(DataflowProcessor processor, List<Object> MOCK_MESSAGES) {
            log.trace "After run > ${name}"
            finalizeTask()
        }

        @Override
        public void afterStop(final DataflowProcessor processor) {
            // increment the session sync
            def val = session.taskDeregister()
            log.debug "After stop > deregister phaser '$name' :: $val"
        }

        @Override
        public Object messageArrived(final DataflowProcessor processor, final DataflowReadChannel<Object> channel, final int index, final Object message) {
            log.trace "Message arrived > task: ${name}; channel: $index; value: $message"
            return message
        }


        public Object controlMessageArrived(final DataflowProcessor processor, final DataflowReadChannel<Object> channel, final int index, final Object message) {
            log.trace "Control message arrived > task: ${name}; channel: $index; value: $message"

            // intercepts 'poison pill' message e.g. termination control message
            // in order the launch the task execution on the underlying system
            if( message == PoisonPill.instance)  {

                if( mergeIndex.get()>0 ) {

                    def task = MergeTaskProcessor.this.createTaskRun()
                    try {
                        MergeTaskProcessor.this.mergeTaskRun(task)
                    }
                    catch( Throwable e ) {
                        handleException(e, task)
                    }

                }
                else {
                    log.warn "No data collected by task > $name -- Won't execute it. Something may be wrong in your execution flow"
                }

            }
            return message
        }

        public boolean onException(final DataflowProcessor processor, final Throwable e) {
            handleException(e)
        }

    }


}
