package nextflow.processor
import groovy.transform.InheritConstructors
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.operator.DataflowEventAdapter
import groovyx.gpars.dataflow.operator.DataflowOperator
import groovyx.gpars.dataflow.operator.DataflowProcessor
import nextflow.exception.InvalidExitException
import nextflow.util.CacheHelper
import nextflow.util.FileHelper
/**
 * Defines the parallel tasks execution logic
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@Slf4j
@InheritConstructors
class ParallelTaskProcessor extends TaskProcessor {

    /**
     * Keeps track of the task instance executed by the current thread
     */
    protected final ThreadLocal<TaskRun> currentTask = new ThreadLocal<>()

    @Override
    protected void createOperator() {


        def opInputs = new ArrayList(taskConfig.inputs.channels)
        def opOutputs = new ArrayList(taskConfig.outputs.channels)

        /*
         * check if there are some iterators declaration
         * the list holds the index in the list of all *inputs* for the {@code each} declaration
         */
        def iteratorIndexes = []
        taskConfig.inputs.eachWithIndex { param, index ->
            if( param instanceof EachInParam ) {  iteratorIndexes << index }
        }

        /*
         * When one (or more) {@code each} are declared as input, it is created an extra
         * operator which will receive the inputs from the channel (excepts the values over iterate)
         *
         * The operator will *expand* the received inputs, iterating over the user provided value and
         * forwarding the final values the the second *parallel* processor executing the user specified task
         */
        if( iteratorIndexes ) {
            final size = taskConfig.inputs.size()
            // the script implementing the iterating process
            final forwarder = createForwardWrapper(size, iteratorIndexes)
            // the channel forwarding the data from the *iterator* process to the target task
            final linkingChannels = new ArrayList(size)
            size.times { linkingChannels[it] = new DataflowQueue() }

            // instantiate the iteration process
            def params = [inputs: opInputs, outputs: linkingChannels, maxForks: 1, listeners: [new IteratorProcessInterceptor()]]
            session.allProcessors << (processor = new DataflowOperator(group, params, forwarder).start())

            // set as next inputs the result channels of the iteration process
            opInputs = linkingChannels
        }


        /*
         * create a mock closure to trigger the operator
         */
        final wrapper = createCallbackWrapper( taskConfig.inputs.size(), this.&invokeTask )

        /*
         * create the output
         */
        def maxForks = taskConfig.maxForks ?: session.config.poolSize
        log.debug "Creating operator > $name -- maxForks: $maxForks"

        def params = [inputs: opInputs, outputs: opOutputs, maxForks: maxForks, listeners: [new TaskProcessorInterceptor()] ]
        session.allProcessors << (processor = new DataflowOperator(group, params, wrapper).start())


    }

    /**
     * Implements the closure which *combines* all the iteration
     *
     * @param numOfInputs Number of in/out channel
     * @param indexes The list of indexes which identify the position of iterators in the input channels
     * @return The clousre implementing the iteration/forwarding logic
     */
    protected createForwardWrapper( int numOfInputs, List indexes ) {

        final args = []
        numOfInputs.times { args << "x$it" }

        /*
         * Explaining the following closure:
         *
         * - it has to be evaluated as a string since the number of input must much the number input channel
         *   that is known only at runtime
         *
         * - 'out' holds the list of all input values which need to be forwarded (binded) to output as many
         *   times are the items in the iteration list
         *
         * - the iteration list(s) is (are) passed like in the closure inputs like the other values,
         *   the *indexes* argument defines the which of them are the iteration lists
         *
         * - 'itr' holds the list of all iteration lists
         *
         * - using the groovy method a combination of all values is create (cartesian product)
         *   see http://groovy.codehaus.org/groovy-jdk/java/util/Collection.html#combinations()
         *
         * - the resulting values are replaced in the 'out' array of values and forwarded out
         *
         */

        final str =
            """
            { ${args.join(',')} ->
                def out = [ ${args.join(',')} ]
                def itr = [ ${indexes.collect { 'x'+it }.join(',')} ]
                def cmb = itr.combinations()
                for( entries in cmb ) {
                    def count = 0
                    n.times { i->
                        if( i in indexes ) { out[i] = entries[count++] }
                    }
                    bindAllOutputValues( *out )
                }
            }
            """

        final Binding binding = new Binding( indexes: indexes, n: numOfInputs )
        final result = (Closure)new GroovyShell(binding).evaluate (str)

        return result

    }

    /**
     * Create the {@code TaskDef} data structure and initialize the task execution context
     * with the received input values
     *
     * @param values
     * @return
     */
    final protected TaskRun initTaskRun(List values) {

        final TaskRun task = createTaskRun()

        /*
         * initialize the inputs/outputs for this task instance
         */
        taskConfig.inputs.each { InParam param -> task.setInput(param) }
        taskConfig.outputs.each { OutParam param -> task.setOutput(param) }

        // -- map the inputs to a map and use to delegate closure values interpolation
        def map = new DelegateMap(ownerScript)

        /*
         * initialize the inputs for this task instances
         */
        taskConfig.inputs.eachWithIndex { InParam param, int index ->

            // add the value to the task instance
            task.setInput(param, values.get(index))

            // otherwise put in on the map used to resolve the values evaluating the script
            if( param instanceof ValueInParam || param instanceof EachInParam ) {
                map[ param.name ] = values.get(index)
            }

        }

        /*
         * initialize the task code to be executed
         */
        task.code = this.code.clone() as Closure
        task.code.delegate = map
        task.code.setResolveStrategy(Closure.DELEGATE_FIRST)

        return task
    }


    /**
     * The processor execution body
     *
     * @param processor
     * @param values
     */
    final protected void invokeTask( def args ) {

        // create and initialize the task instance to be executed
        List params = args instanceof List ? args : [args]
        final task = initTaskRun(params)

        // -- call the closure and execute the script
        try {

            currentTask.set(task)
            task.script = task.code.call()?.toString()?.stripIndent()

            // create an hash for the inputs and code
            // Does it exist in the cache?
            // NO --> launch the task
            // YES --> Does it exited OK and exists the expected outputs?
            //          YES --> return the outputs
            //          NO  --> launch the task

            def keys = [session.uniqueId, task.script ]
            // add all the input name-value pairs to the key generator
            task.inputs.each { keys << it.key.name << it.value }

            def hash = CacheHelper.hasher(keys).hash()
            def folder = FileHelper.createWorkFolder(session.workDir, hash)
            def cached = session.cacheable && taskConfig.cacheable && checkCachedOutput(task,folder)
            if( !cached ) {
                log.info "Running task > ${task.name}"

                // run the task
                task.workDirectory = createTaskFolder(folder, hash)
                launchTask( task )

                // save the exit code
                new File(folder, '.exitcode').text = task.exitCode

                // check if terminated successfully
                boolean success = (task.exitCode in taskConfig.validExitCodes)
                if ( !success ) {
                    throw new InvalidExitException("Task '${task.name}' terminated with an invalid exit code: ${task.exitCode}")
                }

                // -- bind output (files)
                collectOutputs(task)
                bindOutputs(task)
            }

        }
        finally {
            lastRunTask = task
            task.status = TaskRun.Status.TERMINATED
        }
    }


    class TaskProcessorInterceptor extends DataflowEventAdapter {

        @Override
        public void afterStart(final DataflowProcessor processor) {
            // increment the session sync
            def val = session.taskRegister()
            log.trace "After start > register phaser '$name' :: $val"
        }

        @Override
        public List<Object> beforeRun(final DataflowProcessor processor, final List<Object> messages) {
            log.trace "Before run > ${currentTask.get()?.name ?: name} -- messages: ${messages}"
            return messages;
        }

        @Override
        public Object messageArrived(final DataflowProcessor processor, final DataflowReadChannel<Object> channel, final int index, final Object message) {
            if( log.isTraceEnabled() ) {
                def channelName = taskConfig.inputs?.names?.get(index)
                def taskName = currentTask.get()?.name ?: name
                log.trace "Message arrived > ${taskName} -- ${channelName} => ${message}"
            }

            return message;
        }

        @Override
        public Object controlMessageArrived(final DataflowProcessor processor, final DataflowReadChannel<Object> channel, final int index, final Object message) {
            if( log.isTraceEnabled() ) {
                def channelName = taskConfig.inputs?.names?.get(index)
                def taskName = currentTask.get()?.name ?: name
                log.trace "Control message arrived > ${taskName} -- ${channelName} => ${message}"
            }

            return message;
        }

        @Override
        void afterRun(DataflowProcessor processor, List<Object> MOCK_MESSAGES) {
            log.trace "After run > ${currentTask.get()?.name ?: name}"
            currentTask.remove()
            finalizeTask()
        }

        @Override
        public void afterStop(final DataflowProcessor processor) {
            // increment the session sync
            def val = session.taskDeregister()
            log.debug "After stop > deregister phaser '$name' :: $val"
        }

        /**
         * Invoked if an exception occurs. Unless overridden by subclasses this implementation returns true to terminate the operator.
         * If any of the listeners returns true, the operator will terminate.
         * Exceptions outside of the operator's body or listeners' messageSentOut() handlers will terminate the operator irrespective of the listeners' votes.
         * When using maxForks, the method may be invoked from threads running the forks.
         * @param processor
         * @param error
         * @return
         */
        public boolean onException(final DataflowProcessor processor, final Throwable error) {
            handleException( error, currentTask.get() )
        }

    }

    /*
     * logger class for the *iterator* processor
     */
    class IteratorProcessInterceptor extends DataflowEventAdapter {

        @Override
        public boolean onException(final DataflowProcessor processor, final Throwable e) {
            log.error "task '$name' > error on internal iteration process", e
            return true;
        }

        @Override
        public Object messageArrived(final DataflowProcessor processor, final DataflowReadChannel<Object> channel, final int index, final Object message) {
            log.trace "task '$name' > message arrived for iterator '${taskConfig.inputs.names[index]}' with value: '$message'"
            return message;
        }

        @Override
        public Object messageSentOut(final DataflowProcessor processor, final DataflowWriteChannel<Object> channel, final int index, final Object message) {
            log.trace "task '$name' > message forwarded for iterator '${taskConfig.inputs.names[index]}' with value: '$message'"
            return message;
        }


        @Override
        public Object controlMessageArrived(final DataflowProcessor processor, final DataflowReadChannel<Object> channel, final int index, final Object message) {
            log.trace "task '$name' > control message arrived for iterator '${taskConfig.inputs.names[index]}'"
            return message;
        }
    }


}
