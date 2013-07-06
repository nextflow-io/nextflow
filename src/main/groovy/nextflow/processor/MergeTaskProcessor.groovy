package nextflow.processor
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
 * Defines the 'merge' operation logic
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@Slf4j
@InheritConstructors
class MergeTaskProcessor extends TaskProcessor {


    /*
     * The merge task is composed by two operator, the first creates a single 'script' to be executed by second one
     */
    @Override
    protected void createOperator() {

        log.debug "Starting merge > ${name}"

        final args = []
        taskConfig.inputs.size().times { args << "__p$it" }

        final str = " { ${args.join(',')} -> callback([ ${args.join(',')} ]) }"
        final binding = new Binding( ['callback': this.&mergeScriptCollector] )
        final wrapper = (Closure)new GroovyShell(binding).evaluate (str)

        mergeHashesList = new LinkedList<>()
        mergeTempFolder = FileHelper.createTempFolder()

        def params = [inputs: new ArrayList(taskConfig.inputs.values()), outputs: new ArrayList(taskConfig.outputs.values()), listeners: [new MergeProcessorInterceptor()] ]
        processor = new DataflowOperator(group, params, wrapper)
        session.allProcessors.add(processor)

        // increment the session sync
        session.sync.countUp()

        // -- start it
        processor.start()

    }

    protected void mergeTaskRun(TaskRun task) {

        try {

            // -- create the unique hash number for this tasks,
            //    collecting the id of all the executed runs
            //    and sorting them
            def hasher = CacheHelper.hasher(session.uniqueId)
            mergeHashesList.sort()
            mergeHashesList.each { Integer entry ->  hasher = CacheHelper.hasher(hasher,entry) }
            def hash = hasher.hash()
            log.trace "Merging task > $name -- hash: $hash"

            def folder = FileHelper.createWorkFolder(hash)
            log.trace "Merging task > $name -- trying cached: $folder"

            def cached = session.cacheable && taskConfig['cacheable'] && checkCachedOutput(task,folder)
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
                new File(folder, '.exitcode').text = task.exitCode

                // -- check if terminated successfully
                boolean success = (task.exitCode in taskConfig.validExitCodes)
                if ( !success ) {
                    throw new InvalidExitException("Task '${task.name}' terminated with an invalid exit code: ${task.exitCode}")
                }

                bindOutputs(task)
            }
        }
        finally {
            task.status = TaskRun.Status.TERMINATED
        }

    }

    private List<Integer> mergeHashesList

    private File mergeTempFolder

    private AtomicInteger mergeIndex = new AtomicInteger()

    private def mergeScript = new StringBuilder()

    protected void mergeScriptCollector( List params ) {
        final currentIndex = mergeIndex.incrementAndGet()
        log.info "Collecting task > ${name} ($currentIndex)"

        // -- map the inputs to a map and use to delegate closure values interpolation
        def inputVars = new DelegateMap(ownerScript)
        taskConfig.inputs?.keySet()?.eachWithIndex { name, index ->
            inputVars[name] = params[index]
        }

        /*
         * initialize the task code to be executed
         */
        Closure scriptClosure = this.code.clone() as Closure
        scriptClosure.delegate = inputVars
        scriptClosure.setResolveStrategy(Closure.DELEGATE_FIRST)

        def commandToRun = normalizeScript(scriptClosure.call()?.toString())

        /*
         * create a unique hash-code for this task run and save it into a list
         * which maintains all the hashes for executions making-up this merge task
         */
        def keys = [commandToRun]
        if( inputVars.containsKey('-') ) {
            keys << inputVars['-']
        }
        keys << 7
        mergeHashesList << CacheHelper.hasher(keys).hash().asInt()

        /*
         * save the script to execute into a separate unique-named file
         */
        def index = currentIndex
        def scriptName = ".merge_command.sh.${index.toString().padLeft(4,'0')}"
        def scriptFile = new File(mergeTempFolder, scriptName)
        scriptFile.text = commandToRun
        scriptFile.setExecutable(true)

        // the command to launch this command
        def scriptCommand = scriptFile.absolutePath

        // check if some input have to be send
        if( inputVars.containsKey('-') ) {
            def inputName = ".merge_command.input.$index"
            def inputFile = new File( mergeTempFolder, inputName )
            inputFile.text = inputVars['-']

            // pipe the user input to the user command
            scriptCommand = "$scriptCommand < ${inputFile.toString()}"
        }

        // create a unique script collecting all the commands
        mergeScript << scriptCommand << '\n'

    }



    /**
     * A task of type 'merge' binds the output when it terminates it's work, i.e. when
     * it receives a 'poison pill message that will stop it
     */
    class MergeProcessorInterceptor extends DataflowEventAdapter {

        @Override
        void afterRun(DataflowProcessor processor, List<Object> MOCK_MESSAGES) {
            log.trace "After run > ${name}"
            finalizeTask()
        }

        @Override
        public void afterStop(final DataflowProcessor processor) {
            log.debug "After stop > $name"
            // increment the session sync
            session.sync.countDown()
        }

        @Override
        public Object messageArrived(final DataflowProcessor processor, final DataflowReadChannel<Object> channel, final int index, final Object message) {
            log.trace "Received message > task: ${name}; channel: $index; value: $message"
            return message
        }

        @Override
        public void afterStart(final DataflowProcessor processor) {
            log.trace "After start > $name"
        }


        public Object controlMessageArrived(final DataflowProcessor processor, final DataflowReadChannel<Object> channel, final int index, final Object message) {

            // intercepts 'poison pill' message e.g. termination control message
            // in order the launch the task execution on the underlying system
            if( message == PoisonPill.instance)  {

                def task = createTaskRun()
                try {
                    mergeTaskRun(task)
                }
                catch( Exception e ) {
                    handleException(e, task)
                }

            }
            return message
        }

        public boolean onException(final DataflowProcessor processor, final Throwable e) {
            handleException(e)
        }


    }


}
