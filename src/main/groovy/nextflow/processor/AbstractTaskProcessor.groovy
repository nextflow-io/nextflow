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
import java.util.concurrent.atomic.AtomicInteger
import java.util.concurrent.locks.ReentrantLock

import com.google.common.hash.HashCode
import groovy.io.FileType
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.Dataflow
import groovyx.gpars.dataflow.DataflowBroadcast
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowVariable
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.operator.DataflowEventAdapter
import groovyx.gpars.dataflow.operator.DataflowOperator
import groovyx.gpars.dataflow.operator.DataflowProcessor
import groovyx.gpars.dataflow.operator.PoisonPill
import groovyx.gpars.dataflow.stream.DataflowStreamWriteAdapter
import groovyx.gpars.group.PGroup
import nextflow.Nextflow
import nextflow.Session
import nextflow.exception.InvalidExitException
import nextflow.exception.MissingFileException
import nextflow.exception.TaskValidationException
import nextflow.script.AbstractScript
import nextflow.util.CacheHelper
import nextflow.util.FileHelper
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
abstract class AbstractTaskProcessor implements TaskProcessor {

    static final File TMP_FOLDER = new File('tmp')

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
    protected int maxForks = 1

    protected PGroup group = Dataflow.retrieveCurrentDFPGroup()

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
     * The system interpreter used to execute the script, by default {@code 'bash'}
     */
    protected shell = ['bash','-ue']

    /**
     * The exit code which define a valid result, default {@code 0}
     */
    protected List<Integer> validExitCodes = [0]

    protected ErrorStrategy errorStrategy = ErrorStrategy.TERMINATE

    /**
     * When {@code true} the task stdout is redirected to the application console
     */
    protected boolean echo

    /**
     * When the task is cacheable, by default {@code true}
     */
    protected boolean cacheable = true


    // ---== private ==---

    protected TaskRun lastRunTask

    protected DataflowProcessor processor

    private allScalarValues

    private static final creationLock = new ReentrantLock(true)

    private static final folderLock = new ReentrantLock(true)

    private final random = new Random()

    private Boolean errorShown = Boolean.FALSE

    private final errorLock = new ReentrantLock(true)

    AbstractTaskProcessor( Session session ) {
        this.session = session
    }

    AbstractTaskProcessor( Session session, AbstractScript script,  boolean bindOnTermination ) {
        this.session = session
        this.ownerScript = script
        this.bindOnTermination = bindOnTermination

        this.maxForks = bindOnTermination ? 1 : session.config.poolSize

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
    AbstractTaskProcessor shell( value ) {
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

    AbstractTaskProcessor errorStrategy( String value ) {
        this.errorStrategy = ErrorStrategy.valueOf(value?.toUpperCase())
        return this
    }


    @Override
    AbstractTaskProcessor cacheable( boolean value ) {
        this.cacheable = value
        return this
    }

    @Override
    AbstractTaskProcessor maxForks(int max) {
        this.maxForks = max
        return this
    }

    @Override
    AbstractTaskProcessor script(Closure script) {
        this.code = script
        return this
    }

    @Override
    AbstractTaskProcessor script( def shell, Closure script) {
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
    int getMaxForks() { return maxForks }

    @Override
    def getShell() { shell }

    @Override
    List<Integer> getValidExitCodes() { validExitCodes }

    @Override
    ErrorStrategy getErrorStrategy() { errorStrategy }

    @Override
    boolean getCacheable() { cacheable }

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

        /*
         * generate the processor name if not specified
         */
        if ( !name ) {
            name = "task${session.allProcessors.size()+1}"
        }

        /*
         * Normalize the input channels:
         * - at least one input channel have to be provided,
         *   if missing create an dummy 'input' set to true
         */

        if( inputs.size() == 0 ) {
            input('$':true)
        }

        allScalarValues = !inputs.values().any { !(it instanceof DataflowVariable) }

        /*
         * Normalize the output
         * - even though the output may be empty, let return the stdout as output by default
         */
        if ( outputs.size() == 0 ) {
            output('-')
        }

        /*
         * bind the outputs to the script scope
         */
        if( ownerScript ) {
            outputs.each { name, channel ->
                if( name != '-' ) { ownerScript.setProperty(name, channel) }
            }
        }

        if( bindOnTermination ) {
            mergeCreateOperator()
        }
        else {
            createTaskOperator()
        }


        /*
         * When there is a single output channel, return let returns that item
         * otherwise return the list
         */
        def result = outputs.values()
        return result.size() == 1 ? result[0] : result
    }


    private void createTaskOperator() {

        def opInputs = new ArrayList(inputs.values())
        def opOutputs = new ArrayList(outputs.values())


        /*
         * create a mock closure to trigger the operator
         */
        Closure mock = createMockClosure()

        /*
         * create the output
         */
        def params = [inputs: opInputs, outputs: opOutputs, maxForks: maxForks, listeners: [new TaskProcessorInterceptor()] ]
        session.allProcessors << (processor = new DataflowOperator(group, params, mock).start())
        // increment the session sync
        session.sync.countUp()

    }

    /*
     * The merge task is composed by two operator, the first creates a single 'script' to be executed by second one
     */
    private void mergeCreateOperator() {

        log.debug "Starting merge > ${name}"

        final args = []
        inputs.size().times { args << "__p$it" }

        final str = " { ${args.join(',')} -> callback([ ${args.join(',')} ]) }"
        final binding = new Binding( ['callback': this.&mergeScriptCollector] )
        final wrapper = (Closure)new GroovyShell(binding).evaluate (str)

        mergeHashesList = new LinkedList<>()
        mergeTempFolder = FileHelper.createTempFolder()

        def params = [inputs: new ArrayList(inputs.values()), outputs: new ArrayList(outputs.values()), listeners: [new MergeProcessorInterceptor()] ]
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

            def cached = session.cacheable && this.cacheable && checkCachedOutput(task,folder)
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
                boolean success = (task.exitCode in validExitCodes)
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

    protected String getShellCommandString() {
        shell instanceof List ? shell.join(' ') : shell?.toString()
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
        inputs?.keySet()?.eachWithIndex { name, index ->
            inputVars[name] = params[index]
        }

        /*
         * initialize the task code to be executed
         */
        Closure scriptClosure = this.code.clone() as Closure
        scriptClosure.delegate = inputVars
        scriptClosure.setResolveStrategy(Closure.DELEGATE_FIRST)

        def commandToRun = scriptClosure.call()?.toString()?.stripIndent()?.trim()

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

        // the command to launch this command
        def scriptCommand =  getShellCommandString() + ' ' + scriptFile

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


    Closure createMockClosure() {

        // the closure by the GPars dataflow constraints MUST have has many parameters
        // are the input channels, so let create it on-fly
        final params = []
        inputs.size().times { params << "__\$$it" }
        final str = "{ ${params.join(',')} -> runTask(__\$0) }"

        final binding = new Binding( ['runTask': this.&runTask] )
        (Closure)new GroovyShell(binding).evaluate (str)
    }


    final protected TaskRun createTaskRun() {

        TaskRun result
        creationLock.lock()
        try {
            def num = session.tasks.size()
            result = new TaskRun(id: num, status: TaskRun.Status.PENDING, index: ++index, name: "$name ($index)" )
            session.tasks.put( this, result )
        }
        finally {
            creationLock.unlock()
        }

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
        log.debug "Creating a new task > $name"

        final TaskRun task = createTaskRun()

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
            def cached = session.cacheable && this.cacheable && checkCachedOutput(task,folder)
            if( !cached ) {
                log.info "Running task > ${task.name}"

                // run the task
                task.workDirectory = createTaskFolder(folder, hash)
                launchTask( task )

                // save the exit code
                new File(folder, '.exitcode').text = task.exitCode

                // check if terminated successfully
                boolean success = (task.exitCode in validExitCodes)
                if ( !success ) {
                    throw new InvalidExitException("Task '${task.name}' terminated with an invalid exit code: ${task.exitCode}")
                }

                // -- bind output (files)
                if ( success && !bindOnTermination ) {
                    bindOutputs(task)
                }
            }

        }
        finally {
            lastRunTask = task
            task.status = TaskRun.Status.TERMINATED
        }
    }

    final protected File createTaskFolder( File folder, HashCode hash  ) {

        folderLock.lock()
        try {

            if( folder.exists() ) {

                // find another folder name that does NOT exist
                while( true ) {
                    hash = CacheHelper.hasher( [hash.asInt(), random.nextInt() ] ).hash()
                    folder = FileHelper.createWorkFolder(hash)
                    if( !folder.exists() ) {
                        break
                    }
                }
            }

            if( !folder.mkdirs() ) {
                throw new IOException("Unable to create folder: $folder -- check file system permission")
            }

            return folder
        }

        finally {
            folderLock.unlock()
        }

    }


    final checkCachedOutput(TaskRun task, File folder) {
        if( !folder.exists() ) {
            log.trace "Cached folder does not exists > $folder -- return false"
            // no folder -> no cached result
            return false
        }

        def exitFile = new File(folder,'.exitcode')
        if( exitFile.isEmpty() ) {
            log.trace "Exit file is empty > $exitFile -- return false"
            return false
        }

        def exitValue = exitFile.text.trim()
        def exitCode = exitValue.isInteger() ? exitValue.toInteger() : null
        if( exitCode == null || !(exitCode in validExitCodes) ) {
            log.trace "Exit code is not valid > $exitValue -- return false"
            return false
        }

        try {
            task.exitCode = exitCode
            task.workDirectory = folder

            // -- check if all output resources are available
            def producedFiles = collectAndValidateOutputs(task)
            log.info "Cached task > ${task.name}"

            // -- print out the cached tasks output when 'echo' is true
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

    /**
     * Handles an error raised during the processor execution
     *
     * @param error The exception raised during the task execution
     * @param task The {@code TaskDef} instance which raised the exception
     * @return {@code true} to terminate the processor execution,
     *         {@code false} ignore the error and continue to process other pending tasks
     */
    final protected boolean handleException( Throwable error, TaskRun task = null ) {

        // -- when is a task level error and the user has chosen to ignore error, just report and error message
        //    return 'false' to DO NOT stop the execution
        if( error instanceof TaskValidationException && errorStrategy == ErrorStrategy.IGNORE ) {
            log.warn "Error running task > ${error.getMessage()} -- error is ignored"
            return false
        }

        // -- synchronize on the errorLock to avoid multiple report of the same error
        errorLock.lock()
        try {
            if( errorShown ) { return true }
            errorShown = Boolean.TRUE
        }
        finally {
            errorLock.unlock()
        }


        if( error instanceof TaskValidationException ) {

            // compose a readable error message
            def message = []
            message << error.getMessage()

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

            message << "\nTip: when you have fixed the problem you may continue the execution appending to the nextflow command line the '-resume' option"

            log.error message.join('\n')

        }
        else {
            log.error("Failed to execute task > '${task?.name ?: name}' -- Check the log file '.nextflow.log' for more details", error)
        }


        session.abort()
        return true
    }

    final protected void finalizeTask( ) {

        final task = currentTask.get()
        currentTask.remove()

        // -- when all input values are 'scalar' (not queue or broadcast stream)
        //    stops after the first run
        if( allScalarValues ) {
            log.debug "Finalize > ${task?.name ?: name} terminates since all values are scalar"
            // send a poison pill to 'termiate' all downstream channels
            sendPoisonPill()
            processor.terminate()
        }

    }

    final protected synchronized void sendPoisonPill() {

        outputs.each { name, channel ->
            if( channel instanceof DataflowQueue ) {
                log.trace "Sending Poison-pill over $name channel"
                channel.bind( PoisonPill.instance )
            }
            else if( channel instanceof DataflowStreamWriteAdapter ) {
                log.trace "Sending Poison-pill over $name channel"
                channel.bind( PoisonPill.instance )
            }
            else {
                log.trace "Poison pill is not sent over $name channel"
            }

        }

    }

    /**
     * Bind the expected output files to the corresponding output channels
     * @param processor
     */
    protected void bindOutputs( TaskRun task ) {

        bindOutputs(collectAndValidateOutputs(task))

    }

    protected Map collectAndValidateOutputs( TaskRun task ) {

        // -- collect the produced output
        def allFiles = [:]
        outputs.keySet().each { name ->
            allFiles[ name ] = collectResultFile(task, name)
        }

        return allFiles
    }

    synchronized protected bindOutputs( Map allOutputResources ) {

        log.trace "Binding results > task: ${currentTask.get()?.name ?: name} - values: ${allOutputResources}"

        // -- bind each produced file to its own channel
        outputs.keySet().eachWithIndex { name, index ->

            def entry = allOutputResources[name]

            if( name == '-' && entry instanceof File ) {
                processor.bindOutput(index, entry.text)
            }
            else if( entry instanceof Collection ) {
                entry.each { processor.bindOutput(index, it) }
            }
            else {
                processor.bindOutput(index, entry)
            }
        }
    }

    private class BaseEventHandler extends DataflowEventAdapter {

        @Override
        void afterRun(DataflowProcessor processor, List<Object> MOCK_MESSAGES) {
            log.trace "After run > ${currentTask.get()?.name ?: name}"
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

    }

    class TaskProcessorInterceptor extends BaseEventHandler {

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


    /**
     * A task of type 'merge' binds the output when it terminates it's work, i.e. when
     * it receives a 'poison pill message that will stop it
     */
    class MergeProcessorInterceptor extends BaseEventHandler {

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
     * Execute the specified task shell script
     *
     * @param script The script string to be execute, e.g. a BASH script
     * @return {@code TaskDef}
     */
    protected abstract void launchTask( TaskRun task )

    /**
     * The file which contains the stdout produced by the executed task script
     *
     * @param task The user task to be executed
     * @return The absolute file to the produced script output
     */
    protected abstract getStdOutFile( TaskRun task )

    /**
     * Collect the file(s) with the name specified, produced by the execution
     *
     * @param path The job working path
     * @param fileName The file name, it may include file name wildcards
     * @return The list of files matching the specified name
     */
    protected collectResultFile( TaskRun task, String fileName ) {
        assert fileName
        assert task
        assert task.workDirectory

        // the '-' stands for the script stdout, save to a file
        if( fileName == '-' ) {
            return getStdOutFile(task)
        }

        // replace any wildcards characters
        // TODO give a try to http://code.google.com/p/wildcard/  -or- http://commons.apache.org/io/
        String filePattern = fileName.replace("?", ".?").replace("*", ".*?")

        // when there's not change in the pattern, try to find a single file
        if( filePattern == fileName ) {
            def result = new File(task.workDirectory,fileName)
            if( !result.exists() ) {
                throw new MissingFileException("Missing output file: '$fileName' expected by task: ${this.name}")
            }
            return result
        }

        // scan to find the file with that name
        List files = []
        task.workDirectory.eachFileMatch(FileType.FILES, ~/$filePattern/ ) { File it -> files << it}
        if( !files ) {
            throw new MissingFileException("Missing output file(s): '$fileName' expected by task: ${this.name}")
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

