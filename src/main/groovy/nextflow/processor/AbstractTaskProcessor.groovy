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

import com.google.common.hash.HashCode
import com.google.common.hash.Hashing
import groovy.io.FileType
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.Dataflow
import groovyx.gpars.dataflow.DataflowBroadcast
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowVariable
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.operator.DataflowEventAdapter
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
import nextflow.util.CacheHelper
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

    private Boolean errorShown = Boolean.FALSE

    /**
     * When {@code true} the task stdout is redirected to the application console
     */
    protected boolean echo

    // ---== private ==---

    protected TaskDef lastRunTask

    protected DataflowProcessor processor

    private allScalarValues

    private creationLock = new ReentrantLock(true)


    private static final folderLock = new ReentrantLock(true)

    private File sharedFolder

    private random = new Random()

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
                this.inputs.put( name, Nextflow.channel(value) )
            }
            // wrap any array with a DataflowQueue
            else if ( value && value.class.isArray() ) {
                this.inputs.put( name, Nextflow.channel(value as List) )
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
            files.each { name -> outputs.put( name, Nextflow.channel() ) }
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
        def params = [inputs: _inputs, outputs: _outputs, maxForks: threads, listeners: [new DataflowInterceptor()] ]
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

        // the closure by the GPars dataflow constraints MUST have has many parameters
        // are the input channels, so let create it on-fly
        final params = []
        inputs.size().times { params << "__\$$it" }
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
    final protected TaskDef initTaskRun(List values) {
        log.debug "Creating a new task > $name"

        final TaskDef task = null
        creationLock.lock()
        try {
            def num = session.tasks.size()
            task = new TaskDef(id: num, status: TaskDef.Status.PENDING, index: ++index, name: "$name ($index)" )
            session.tasks.put( this, task )
        }
        finally {
            creationLock.unlock()
        }

        // -- map the inputs to a map and use to delegate closure values interpolation
        def map = new DelegateMap(ownerScript)

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


    protected final ThreadLocal<TaskDef> currentTask = new ThreadLocal<>()

    /**
     * The processor execution body
     *
     * @param processor
     * @param values
     */
    final protected runTask(TaskDef task) {
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

            def hash = CacheHelper.hasher( [task.script, task.input, task.code.delegate] ).hash()
            def folder = shareWorkDir && sharedFolder ? sharedFolder : folderForHash(hash)
            def cached = session.cacheable && (!shareWorkDir) && checkCachedOutput(task,folder)
            if( !cached ) {
                log.info "Running task > ${task.name}"
                log.debug "New key: ${hash.toString()}\ntask: ${task.name}\nfolder: ${folder}\nshared: ${sharedFolder}\nscript: ${task.script}\ninput: ${task.input}\nmap: ${task.code.delegate?.dump()}"

                task.workDirectory = createTaskFolder(folder, hash)
                launchTask( task )
            }
            else {
                log.debug "Cached key: ${hash.toString()}\ntask: ${task.name}\nfolder: ${folder}\nshared: ${sharedFolder}\nscript: ${task.script}\ninput: ${task.input}\nmap: ${task.code.delegate?.dump()}"
            }

            boolean success = (task.exitCode in validExitCodes)
            if ( !success ) {
                throw new InvalidExitException("Task '${task.name}' terminated with an invalid exit code: ${task.exitCode}")
            }

            // -- bind output (files)
            if ( success && !bindOnTermination ) {
                bindOutputs(task)
            }
        }
        finally {
            lastRunTask = task
            task.status = TaskDef.Status.TERMINATED
        }
    }


    final protected File folderForHash(HashCode hash) {

        def bucket = Hashing.consistentHash(hash, 100)
        def folder = new File("./tmp/${session.scriptName}/${bucket}", hash.toString()).absoluteFile

        return folder
    }

    final protected File createTaskFolder( File folder, HashCode hash ) {

        folderLock.lock()
        try {
            if( shareWorkDir && folder.exists() ) {
                return folder
            }

            if( folder.exists() ) {

                // find another folder name that does NOT exist
                while( true ) {
                    hash = CacheHelper.hasher( [hash.asInt(), random.nextInt() ] ).hash()
                    folder = folderForHash(hash)
                    if( !folder.exists() ) {
                        break
                    }
                }
            }

            if( !folder.mkdirs() ) {
                throw new IOException("Unable to create folder: $folder -- check file system permission")
            }

            if( shareWorkDir ) sharedFolder = folder

            return folder
        }

        finally {
            folderLock.unlock()
        }

    }

    final checkCachedOutput(TaskDef task, File folder) {
        log.debug "Checking cache folder > task: ${task.name} -- folder: ${folder}"
        if( !folder.exists() ) {
            log.debug "Cached folder does not exists > $folder"
            // no folder -> no cached result
            return false
        }

        def exitFile = new File(folder,'.exitcode')
        if( exitFile.isEmpty() ) {
            log.debug "Cached exit does not exists > $folder"
            return false
        }

        def exitValue = exitFile.text.trim()
        if( !exitValue.isInteger() || !(exitValue.toInteger() in validExitCodes) ) {
            log.debug "Cached invalid exit > $folder"
            return false
        }

        try {
            task.exitCode = exitValue.toInteger()
            task.workDirectory = folder

            // -- check if all output resources are available
            def producedFiles = collectOutputs(task)
            log.info "Cached task > ${task.name} -- folder: $folder"

            // -- print out the program output is required
            if( producedFiles.containsKey('-') && echo ) {
                def out = producedFiles['-']
                if( out instanceof File )  {
                    System.out.print(out.text)
                }
                else if( out ) {
                    System.out.print(out.toString())
                }
            }

            // -- now bind the results
            bindOutputs(producedFiles)

            return true
        }
        catch( MissingFileException e ) {
            log.debug "Missing cached file > ${e.getMessage()} -- folder: $folder"
            task.exitCode = Integer.MAX_VALUE
            task.workDirectory = null
            return false
        }


    }

    final synchronized protected handleException( Throwable e, TaskDef task = null ) {

        // -- synchronize on the errorShow to avoid multiple report of the same task
        if( errorShown ) { return }
        errorShown = Boolean.TRUE

        if( e instanceof TaskValidationException ) {
            if ( errorStrategy == ErrorStrategy.IGNORE ) {
                log.debug "${e.getMessage()} -- error is ignored"
                return
            }

            // compose a readable error message
            def message = []
            message << e.getMessage()

            if( task ) {
                // print the executed command
                message << "Command executed:"
                task.script.eachLine {
                    message << "  $it"
                }

                // if the echo was disabled show, the program out when there's an error
                message << "\nCommand output:"
                task.output.eachLine {
                    message << "  $it"
                }

                message << "\nCommand work dir:\n  ${task.workDirectory}"
            }

            log.error message.join('\n')

        }
        else {
            log.error("Failed to execute task > '${task?.name ?: name}' -- Check the log file '.nextflow.log' for more details", e)
        }


        session.abort()

    }

    final protected void finalizeTask( ) {

        final task = currentTask.get()
        currentTask.remove()

        // -- when all input values are 'scalar' (not queue or broadcast stream)
        //    stops after the first run
        if( allScalarValues ) {
            log.debug "Finalize > ${task?.name ?: name} terminates since all values are scalar"
            processor.terminate()
        }

    }

    /**
     * Bind the expected output files to the corresponding output channels
     * @param processor
     */
    protected synchronized void bindOutputs( TaskDef task ) {

        bindOutputs(collectOutputs(task))

    }

    protected Map collectOutputs( TaskDef task ) {
        // -- collect the produced output
        def allFiles = [:]
        outputs.keySet().each { name ->
            if( name == '-' ) {
                def result = collectResultFile(task, name)
                allFiles[ name ] = result instanceof File ? result.text : result?.toString()
            }
            else {
                allFiles[ name ] = collectResultFile(task, name)
            }
        }

        return allFiles
    }

    protected bindOutputs( Map allFiles ) {
        // -- bind each produced file to its own channel
        outputs.keySet().eachWithIndex { fileName, index ->

            def entry = allFiles[fileName]
            if( entry instanceof Collection ) {
                entry.each { processor.bindOutput(index, it) }
            }
            else {
                processor.bindOutput(index, entry)
            }
        }
    }


    class DataflowInterceptor extends DataflowEventAdapter {

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
            log.trace "After run > ${currentTask.get()?.name ?: name}"
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
            log.trace "Received control message > task: ${currentTask.get()?.name ?: name}; channel: $index; value: $message"

            if ( message == PoisonPill.instance && AbstractTaskProcessor.this.bindOnTermination && AbstractTaskProcessor.this.lastRunTask ) {
                log.debug "Bind on termination > task: ${currentTask.get()?.name ?: name}"

                try {
                    bindOutputs(lastRunTask)
                }
                catch( Throwable error ) {
                    handleException(error)
                }
                lastRunTask = null
            }

            return message
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
            log.trace "Received message > task: ${currentTask.get()?.name ?: name}; channel: $index; value: $message"
            return message
        }

        /**
         * Invoked immediately after the operator starts by a pooled thread before the first message is obtained
         */
        @Override
        public void afterStart(final DataflowProcessor processor) {
            log.trace "After start > $name"
        }

        /**
         * Invoked immediately after the operator terminates
         *
         * @param processor The reporting dataflow operator/selector
         */
        @Override
        public void afterStop(final DataflowProcessor processor) {
            log.trace "After stop > $name"
            // increment the session sync
            session.sync.countDown()
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
    protected abstract void launchTask( TaskDef task )

    /**
     * Collect the file(s) with the name specified, produced by the execution
     *
     * @param path The job working path
     * @param name The file name, it may include file name wildcards
     * @return The list of files matching the specified name
     */

    protected collectResultFile( TaskDef task, String name ) {
        assert name
        assert task
        assert task.workDirectory

        // replace any wildcards characters
        // TODO give a try to http://code.google.com/p/wildcard/  -or- http://commons.apache.org/io/
        String filePattern = name.replace("?", ".?").replace("*", ".*?")

        if( filePattern == name ) {
            def result = new File(task.workDirectory,name)
            if( !result.exists() ) {
                throw new MissingFileException("Missing output file: '$name' expected by task: ${this.name}")
            }
            return result
        }

        // scan to find the file with that name
        List files = []
        task.workDirectory.eachFileMatch(FileType.FILES, ~/$filePattern/ ) { File it -> files << it}
        if( !files ) {
            throw new MissingFileException("Missing output file(s): '$name' expected by task: ${this.name}")
        }

        return files
    }




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
                    log.trace "Unable to find a value for: '\$${property}' on script context"
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

