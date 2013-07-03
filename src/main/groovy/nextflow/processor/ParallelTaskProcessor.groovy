package nextflow.processor

import groovy.transform.InheritConstructors
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowReadChannel
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


    @Override
    protected void createOperator() {

        def opInputs = new ArrayList(taskConfig.inputs.values())
        def opOutputs = new ArrayList(taskConfig.outputs.values())


        /*
         * create a mock closure to trigger the operator
         */
        Closure mock = createMockClosure()

        /*
         * create the output
         */
        def maxForks = taskConfig['maxForks'] ?: session.config.poolSize
        def params = [inputs: opInputs, outputs: opOutputs, maxForks: maxForks, listeners: [new TaskProcessorInterceptor()] ]
        session.allProcessors << (processor = new DataflowOperator(group, params, mock).start())

        // increment the session sync
        session.sync.countUp()

    }



    Closure createMockClosure() {

        // the closure by the GPars dataflow constraints MUST have has many parameters
        // are the input channels, so let create it on-fly
        final params = []
        taskConfig.inputs.size().times { params << "__\$$it" }
        final str = "{ ${params.join(',')} -> runTask(__\$0) }"

        final binding = new Binding( ['runTask': this.&runTask] )
        (Closure)new GroovyShell(binding).evaluate (str)
    }

    /**
     * Create the {@code TaskDef} data structure and initialize the task execution context
     * with the received input values
     *
     * @param values
     * @return
     */
    final protected TaskRun initTaskRun(List values) {
        log.debug "Creating a new task > $name"

        final TaskRun task = createTaskRun()

        // -- map the inputs to a map and use to delegate closure values interpolation
        def map = new DelegateMap(ownerScript)

        taskConfig.inputs?.keySet()?.eachWithIndex { String name, int index ->

            // when the name define for this 'input' is '-'
            // copy the value to the task 'input' attribute
            // it will be used to pipe it to the process stdin
            if( name == '-' ) {
                task.input = values.get(index)
            }

            // otherwise put in on the map used to resolve the values evaluating the script
            else {
                map[ name ] = values.get(index)
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


    protected final ThreadLocal<TaskRun> currentTask = new ThreadLocal<>()

    /**
     * The processor execution body
     *
     * @param processor
     * @param values
     */
    final protected runTask(TaskRun task) {
        assert task

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

            def hash = CacheHelper.hasher( [session.uniqueId, task.script, task.input, task.code.delegate] ).hash()
            def folder = FileHelper.createWorkFolder(hash)
            def cached = session.cacheable && taskConfig['cacheable'] && checkCachedOutput(task,folder)
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
                if ( success ) {
                    bindOutputs(task)
                }
            }

        }
        finally {
            lastRunTask = task
            task.status = TaskRun.Status.TERMINATED
        }
    }


    class TaskProcessorInterceptor extends DataflowEventAdapter {

        @Override
        void afterRun(DataflowProcessor processor, List<Object> MOCK_MESSAGES) {
            log.trace "After run > ${currentTask.get()?.name ?: name}"
            currentTask.remove()
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
            log.trace "Received message > task: ${currentTask.get()?.name ?: name}; channel: $index; value: $message"
            return message
        }

        @Override
        public void afterStart(final DataflowProcessor processor) {
            log.trace "After start > $name"
        }


        /**
         * Invoked when all messages required to trigger the operator become available in the input channels.
         *
         * @param processor
         * @param messages
         * @return
         */
        @Override
        List<Object> beforeRun(DataflowProcessor processor, List<Object> messages) {

            // - prepare and initialize the data structure before execute the task
            // - set the current task parameter on a Thread local variable
            final task = initTaskRun(messages)
            log.trace "Before run > ${task.name} -- messages: $messages"

            // HERE COMES THE HACK !
            // the result 'messages' is used to pass the task instance in the 'mock' closure
            // since there MUST be at least one input parameters it is sure to have at least one element.
            // This is used to pass the TASK instanced as argument, others are ignores
            // See 'createMockClosure'
            messages = new Object[ messages.size() ]
            messages[0] = task

            return messages
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


}
