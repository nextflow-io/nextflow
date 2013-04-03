/*
 * Copyright (c) 2012, the authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
 */
package nextflow.processor
import java.util.concurrent.locks.ReentrantLock

import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.Dataflow
import groovyx.gpars.dataflow.DataflowBroadcast
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowVariable
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.operator.DataflowEventAdapter
import groovyx.gpars.dataflow.operator.DataflowEventListener
import groovyx.gpars.dataflow.operator.DataflowOperator
import groovyx.gpars.dataflow.operator.DataflowProcessor
import groovyx.gpars.dataflow.operator.PoisonPill
import groovyx.gpars.group.PGroup
import nextflow.Nextflow
import nextflow.Session
import nextflow.exception.InvalidExitException
import nextflow.exception.MissingFileException
import nextflow.exception.TaskValidationException
import nextflow.script.AbstractScript
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
abstract class AbstractTaskProcessor implements TaskProcessor {

    protected int index

    /**
     * The current workflow execution session
     */
    protected final Session session

    /**
     * The script object which defines this task
     */
    protected final AbstractScript ownerScript

    /**
     * When {@code true} the output is produced only at the process termination
     */
    protected final bindOnTermination

    /**
     * Task specific environment variable to be used
     */
    protected Map<String,String> environment

    /**
     * All inputs channels for this task definition
     */
    protected Map<String,DataflowReadChannel> inputs = new LinkedHashMap<>()

    /**
     * All outputs channels for this task definition
     */
    protected Map<String,DataflowWriteChannel> outputs = new LinkedHashMap<>()

    /**
     * Maximum number of thread that can be used by this processor
     */
    protected int threads = 1

    protected PGroup group = Dataflow.DATA_FLOW_GROUP

    /**
     * The processor descriptive name
     */
    protected String name

    /**
     * When {@code true} all the tasks launched by the same processor will share the same work directory
     */
    protected Boolean shareWorkDir = Boolean.FALSE

    /**
     * The closure wrapping the script to be executed
     */
    protected Closure code

    /**
     * The system interpreter used to execute the script, by default {@code bash}
     */
    protected String shell = 'bash'

    /**
     * The exit code which define a valid result, default {@code 0}
     */
    protected List<Integer> validExitCodes = [0]


    protected ErrorStrategy errorStrategy = ErrorStrategy.TERMINATE

    /**
     * When {@code true} the task stdout is redirected to the application console
     */
    protected boolean echo

    // ---== private ==---

    protected TaskDef lastRunTask

    private allScalarValues

    private creationLock = new ReentrantLock(true)

    private DataflowProcessor processor

    AbstractTaskProcessor( Session session ) {
        this.session = session
    }

    AbstractTaskProcessor( Session session, AbstractScript script,  boolean bindOnTermination ) {
        this.session = session
        this.ownerScript = script
        this.bindOnTermination = bindOnTermination

        // by definition when 'bindOnTermination' is true all tasks must share
        // the same working directory
        if ( bindOnTermination ) {
            this.shareWorkDir = true
        }
    }

    @Override
    AbstractTaskProcessor environment(Map<String, String> environment) {
        assert environment
        this.environment = new HashMap<>(environment)
        return this
    }

    @Override
    AbstractTaskProcessor input(Map<String,?> inputs) {
        // wrap by a linked map o guarantee the insertion order

        inputs?.each { name, value ->
            if ( value instanceof DataflowBroadcast )  {
                this.inputs.put( name, value.createReadChannel() )
            }
            else if( value instanceof DataflowReadChannel ) {
                this.inputs.put( name, value )
            }
            // wrap any collections with a DataflowQueue
            else if( value instanceof Collection ) {
                this.inputs.put( name, Nextflow.queue(value) )
            }
            // wrap any array with a DataflowQueue
            else if ( value && value.class.isArray() ) {
                this.inputs.put( name, Nextflow.queue(value as List) )
            }
            // wrap a single value with a DataflowVariable
            else {
                this.inputs.put( name, Nextflow.val(value) )
            }
        }


        return this
    }

    @Override
    AbstractTaskProcessor output(String... files) {
        if ( files ) {
            files.each { name -> outputs.put( name, Nextflow.queue() ) }
        }

        return this
    }

    @Override
    AbstractTaskProcessor output(Map<String,DataflowWriteChannel> outputs) {
        if ( outputs ) {
            this.outputs.putAll(outputs)
        }

        return this
    }

    @Override
    AbstractTaskProcessor name( String name ) {
        this.name = name
        return this
    }

    @Override
    AbstractTaskProcessor echo( boolean value ) {
        this.echo = value
        return this
    }

    @Override
    AbstractTaskProcessor shareWorkDir(boolean value) {
        this.shareWorkDir = value
        return this
    }

    @Override
    AbstractTaskProcessor shell( String value ) {
        this.shell = value
        return this
    }

    @Override
    AbstractTaskProcessor validExitCodes( List<Integer> values ) {
        this.validExitCodes = values
        return this
    }

    @Override
    AbstractTaskProcessor errorStrategy( ErrorStrategy value ) {
        this.errorStrategy = value
        return this
    }

    @Override
    AbstractTaskProcessor threads(int max) {
        this.threads = max
        return this
    }

    @Override
    AbstractTaskProcessor script(Closure script) {
        this.code = script
        return this
    }

    @Override
    AbstractTaskProcessor script(String shell, Closure script) {
        this.shell = shell
        this.code = script
        return this
    }


    @Override
    Session getSession() { session }

    boolean getBindOnTermination() {  bindOnTermination  }

    @Override
    String getName() { name }

    @Override
    DataflowReadChannel getInput( String name ) { inputs.get(name) }

    @Override
    DataflowWriteChannel getOutput( String name ) { outputs.get(name) }

    @Override
    boolean getEcho() { return echo }

    @Override
    Map<String,String > getEnvironment() { return environment }

    @Override
    int getThreads() { return threads }

    @Override
    boolean getShareWorkDir() { shareWorkDir }

    @Override
    String getShell() { shell }

    @Override
    List<Integer> getValidExitCodes() { validExitCodes }

    @Override
    ErrorStrategy getErrorStrategy() { errorStrategy }

    /**
     * Launch the 'script' define by the code closure as a local bash script
     *
     * @param code A {@code Closure} retuning a bash script e.g.
     *          <pre>
     *              {
     *                 """
     *                 #!/bin/bash
     *                 do this ${x}
     *                 do that ${y}
     *                 :
     *                 """
     *              }
     *
     * @return {@code this} instance
     */
    @Override
    def run() {

        if ( !code ) {
            throw new IllegalArgumentException("Missing 'script' attribute")
        }

        if( shareWorkDir && threads > 1 ) {
            throw new IllegalArgumentException("The 'shareWorkDirectory' attribute cannot be set TRUE when the specified 'thread's are more than 1")
        }

        /*
         * generate the processor name if not specified
         */
        if ( !name ) {
            name = "task${session.allProcessors.size()}"
        }

        /*
         * Normalize the input channels:
         * - at least one input channel have to be provided,
         *   if missing create an dummy 'input' set to true
         */

        if( inputs.size() == 0 ) {
            input('$':true)
        }
        def _inputs = new ArrayList(inputs.values())

        allScalarValues = !_inputs.any { !(it instanceof DataflowVariable) }

        /*
         * Normalize the output
         * - event though the output may be empty, let return the stdout as output by default
         */
        if ( outputs.size() == 0 ) { output('-') }
        def _outputs = new ArrayList(outputs.values())

        // bind the outputs to the script scope
        if( ownerScript ) {
            outputs.each { name, channel ->
                if( name != '-' ) { ownerScript.setProperty(name, channel) }
            }
        }

        /*
         * create a mock closure to trigger the operator
         */
        Closure mock = createMockClosure()

        /*
         * create the output
         */
        def params = [inputs: _inputs, outputs: _outputs, maxForks: threads, listeners: [createListener()] ]
        session.allProcessors << (processor = new DataflowOperator(group, params, mock).start())
        // increment the session sync
        session.sync.countUp()

        /*
         * When there is a single output channel, return let returns that item
         * otherwise return the list
         */
        return _outputs.size() == 1 ? _outputs[0] : _outputs
    }



    Closure createMockClosure() {

        def str
        // when no input is provided just an empty closure

        // the closure by the GPars dataflow constraints MUST have has many parameters
        // are the input channels, so let create it on-fly
        def params = []
        inputs.size().times { params << "__\$$it" }
        str = "{ ${params.join(',')} -> runTask(__\$0) }"
        log.debug "Mock closure: ${str}"

        def binding = new Binding( ['runTask': this.&runTask] )
        (Closure)new GroovyShell(binding).evaluate (str)
    }


    /**
     * Create the {@code TaskDef} data structure and initialize the task execution context
     * with the received input values
     *
     * @param values
     * @return
     */
    final protected TaskDef initTaskRun(List values) {

        final task = null
        creationLock.lock()
        try {
            def num = session.tasks.size()
            task = new TaskDef(id: num, status: TaskDef.Status.PENDING, index: ++index )
            session.tasks.put( this, task )
        }
        finally {
            creationLock.unlock()
        }

        // -- map the inputs to a map and use to delegate closure values interpolation
        def map = new DelegateMap(ownerScript)
        map['taskIndex'] = task.index
        map['taskId'] = task.id

        inputs?.keySet()?.eachWithIndex { name, index ->

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

    ThreadLocal<TaskDef> currentTask = new ThreadLocal<>()

    /**
     * The processor execution body
     *
     * @param processor
     * @param values
     */
    final protected runTask(TaskDef task) {
        log.debug "Running task: '$name' -- ${task.dump()} "

        // -- call the closure and execute the script
        try {
            currentTask.set(task)
            runScript( task.code.call(), task )

            boolean success = (task.exitCode in validExitCodes)
            if ( !success ) {
                throw new InvalidExitException("Task '$name' terminated with an invalid exit-code: ${task.exitCode}")
            }

            // -- bind output (files)
            if ( success && !bindOnTermination ) {
                bindOutputs(processor, task)
            }
        }
        finally {
            lastRunTask = task
            task.status = TaskDef.Status.TERMINATED
        }

    }

    final protected handleException( Throwable e, TaskDef task = null ) {

        if( e instanceof TaskValidationException ) {
            if ( errorStrategy == ErrorStrategy.IGNORE ) {
                log.debug "${e.getMessage()} -- error is ignored"
                return
            }

            log.error("${e.getMessage()} -- Check the log file '.nextflow.log' for more details", e)

            // if the echo was disabled show, the program out when there's an error
            if ( !echo && task ) {
                println task.output
            }

        }
        else {
            log.error("Failed to execute task: '${name}' -- Check the log file '.nextflow.log' for more details", e)
        }


        session.abort()

    }

    final protected void finalizeTask( ) {

        currentTask.remove()

        // -- when all input values are 'scalar' (not queue or broadcast stream)
        //    stops after the first run
        if( allScalarValues ) {
            log.debug "Processor: '$name' terminates since all values are scalar"
            processor.terminate()
        }

    }

    /**
     * Bind the expected output files to the corresponding output channels
     * @param processor
     */
    protected synchronized void bindOutputs( DataflowProcessor processor, TaskDef task ) {

        // -- collect the produced output
        outputs.keySet().eachWithIndex { fileName, index ->

            if ( fileName == '-' ) {
                processor.bindOutput( index, task.output )
            }

            else {
                def files = collectResultFile( task.workDirectory, fileName )
                if ( !files )  {
                    throw new MissingFileException("Missing output file(s): '$fileName' expected by task: '${this.name}'")
                }

                files .each { File file ->
                    if( file.exists() ) {
                        processor.bindOutput(index, file)
                    }
                    else {
                        throw new MissingFileException("Missing output file(s): '$fileName' expected by task: '${this.name}'")
                    }
                }
            }
        }
    }



    /**
     * The operator listener, it implements some important behavior
     */
    private DataflowEventListener createListener() {

        new DataflowEventAdapter() {

            /**
             * Invoked when all messages required to trigger the operator become available in the input channels.
             *
             * @param processor
             * @param messages
             * @return
             */
            @Override
            List<Object> beforeRun(DataflowProcessor processor, List<Object> messages) {
                if( log.isDebugEnabled() ) {
                    log.debug "Task: '$name' before run"
                }
                else if ( log.isTraceEnabled() ) {
                    log.trace "Task: '$name' before run; messages: $messages"
                }

                // - prepare and initialize the data structure before execute the task
                // - set the current task parameter on a Thread local variable
                final TASK = initTaskRun(messages)

                // HERE COMES THE HACK !
                // the result 'messages' is used to pass the task instance in the 'mock' closure
                // since there MUST be at least one input parameters it is sure to have at least one element.
                // This is used to pass the TASK instanced as argument, others are ignores
                // See 'createMockClosure'
                messages = new Object[ messages.size() ]
                messages[0] = TASK

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
                log.debug "Task: '$name' raised and exception"
                handleException( error, AbstractTaskProcessor.this.currentTask.get() )
                return true
            }

            /**
             * Invoked when the operator completes a single run.
             *
             * @param processor
             * @param MOCK_MESSAGES
             */
            @Override
            void afterRun(DataflowProcessor processor, List<Object> MOCK_MESSAGES) {
                log.debug "Task: '$name' after run"
                finalizeTask()
            }


            /**
             * Invoked when a control message (instances of ControlMessage) becomes available in an input channel.
             *
             * @param processor
             * @param channel
             * @param index
             * @param message
             * @return
             */
            @Override
            public Object controlMessageArrived(final DataflowProcessor processor, final DataflowReadChannel<Object> channel, final int index, final Object message) {
                if( log.isDebugEnabled() ) {
                    log.debug "Task: '$name' control message arrived"
                }
                else if ( log.isTraceEnabled() ) {
                    log.trace "Task: '$name' control message arrived; channel index: $index; value: $message"
                }

                if ( message == PoisonPill.instance && AbstractTaskProcessor.this.bindOnTermination && AbstractTaskProcessor.this.lastRunTask ) {
                    log.debug "Processor '$name' binding on termination"
                    try {
                        bindOutputs(processor, lastRunTask)
                    }
                    catch( Throwable error ) {
                        handleException(error)
                    }
                    lastRunTask = null
                }

                return message;
            }

            /**
             * Invoked when a message becomes available in an input channel.
             *
             * @param processor
             * @param channel
             * @param index
             * @param message
             * @return
             */
            @Override
            public Object messageArrived(final DataflowProcessor processor, final DataflowReadChannel<Object> channel, final int index, final Object message) {
                if( log.isDebugEnabled() ) {
                    log.debug "Task: '$name' receiving message"
                }
                else if ( log.isTraceEnabled() ) {
                    log.trace "Task: '$name' receiving message; channel index: $index; value: $message"
                }

                return message;
            }

            /**
             * Invoked immediately after the operator starts by a pooled thread before the first message is obtained
             */
            @Override
            public void afterStart(final DataflowProcessor processor) {
                log.debug "Task: '$name' after start"
            }

            /**
             * Invoked immediately after the operator terminates
             *
             * @param processor The reporting dataflow operator/selector
             */
            @Override
            public void afterStop(final DataflowProcessor processor) {
                log.debug "Task: '$name' after stop"
                // increment the session sync
                session.sync.countDown()
            }

        }
    }


    protected Map<String,String> getProcessEnvironment() {

        def result = [:]

        // add the config environment entries
        if( session.config.env instanceof Map ) {
            session.config.env.each { name, value ->
                result.put( name, value?.toString() )
            }
        }
        else {
            log.debug "Invalid 'session.config.env' object: ${session.config.env?.class?.name}"
        }

        // add the task specific entries
        if( environment ) {
            environment.each { name, value -> result.put(name, value?.toString()) }
        }

        return result
    }


    /**
     * Execute the specified system script
     *
     * @param script The script string to be execute, e.g. a BASH script
     * @return {@code TaskDef}
     */
    protected abstract void runScript( def script, TaskDef task )

    /**
     * Collect the file(s) with the name specified, produced by the execution
     *
     * @param path The job working path
     * @param name The file name, it may include file name wildcards
     * @return The list of files matching the specified name
     */
    protected abstract List<File> collectResultFile( File path, String name )


    /**
     * Map used to delegate variable resolution to script scope
     */
    @Slf4j
    static class DelegateMap implements Map {

        @Delegate
        private Map<String,Object> local

        private AbstractScript script

        DelegateMap(AbstractScript script) {
            this.script = script
            this.local = [:]
        }

        DelegateMap(Map target) {
            assert target != null
            this.script = script
            this.local = target
        }

        @Override
        public Object get(Object property) {

            if( local.containsKey(property) ) {
                return local.get(property)
            }
            else if ( script ){
                try {
                    return script.getProperty(property?.toString())
                }
                catch( MissingPropertyException e ) {
                    log.debug "Unable to find a value for: '\$${property}' on script context"
                }
            }

            // return the variable name prefixed with the '$' char
            // so give a chance to the bash interpreted to evaluate it
            return '$' + property

        }

        @Override
        public put(String property, Object newValue) {
            local.put(property, newValue)
        }
    }

}

