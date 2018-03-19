/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
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

import static nextflow.processor.ErrorStrategy.*

import java.nio.file.LinkOption
import java.nio.file.NoSuchFileException
import java.nio.file.Path
import java.nio.file.Paths
import java.util.concurrent.atomic.AtomicBoolean
import java.util.concurrent.atomic.AtomicInteger
import java.util.concurrent.atomic.AtomicIntegerArray
import java.util.concurrent.locks.Lock
import java.util.concurrent.locks.ReentrantLock
import java.util.regex.Matcher
import java.util.regex.Pattern

import ch.grengine.Grengine
import com.google.common.hash.HashCode
import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import groovyx.gpars.agent.Agent
import groovyx.gpars.dataflow.Dataflow
import groovyx.gpars.dataflow.DataflowChannel
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.expression.DataflowExpression
import groovyx.gpars.dataflow.operator.DataflowEventAdapter
import groovyx.gpars.dataflow.operator.DataflowOperator
import groovyx.gpars.dataflow.operator.DataflowProcessor
import groovyx.gpars.dataflow.operator.PoisonPill
import groovyx.gpars.dataflow.stream.DataflowStreamWriteAdapter
import groovyx.gpars.group.PGroup
import nextflow.Nextflow
import nextflow.Session
import nextflow.cloud.CloudSpotTerminationException
import nextflow.dag.NodeMarker
import nextflow.exception.FailedGuardException
import nextflow.exception.MissingFileException
import nextflow.exception.MissingValueException
import nextflow.exception.ProcessException
import nextflow.exception.ProcessFailedException
import nextflow.exception.ProcessStageException
import nextflow.exception.ProcessUnrecoverableException
import nextflow.executor.CachedTaskHandler
import nextflow.executor.Executor
import nextflow.extension.DataflowHelper
import nextflow.file.FileHelper
import nextflow.file.FileHolder
import nextflow.file.FilePatternSplitter
import nextflow.script.BaseScript
import nextflow.script.BasicMode
import nextflow.script.EachInParam
import nextflow.script.EnvInParam
import nextflow.script.FileInParam
import nextflow.script.FileOutParam
import nextflow.script.InParam
import nextflow.script.MissingParam
import nextflow.script.OptionalParam
import nextflow.script.OutParam
import nextflow.script.ScriptType
import nextflow.script.SetInParam
import nextflow.script.SetOutParam
import nextflow.script.StdInParam
import nextflow.script.StdOutParam
import nextflow.script.TaskBody
import nextflow.script.ValueInParam
import nextflow.script.ValueOutParam
import nextflow.util.ArrayBag
import nextflow.util.BlankSeparatedList
import nextflow.util.CacheHelper
import nextflow.util.CollectionHelper
/**
 * Implement nextflow process execution logic
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class TaskProcessor {

    /**
     * Implements the closure which *combines* all the iteration
     *
     * @param numOfInputs Number of in/out channel
     * @param indexes The list of indexes which identify the position of iterators in the input channels
     * @return The closure implementing the iteration/forwarding logic
     */
    static class ForwardClosure extends Closure {

        final private Integer len

        final private int numOfParams

        final private List<Integer> indexes

        ForwardClosure(int len, List<Integer> indexes) {
            super(null, null);
            assert len>=1
            this.len = len
            this.numOfParams = len+1
            this.indexes = indexes
        }

        @Override
        public int getMaximumNumberOfParameters() {
            numOfParams
        }

        @Override
        public Class[] getParameterTypes() {
            def result = new Class[numOfParams]
            for( int i=0; i<result.size(); i++ )
                result[i] = Object
            return result
        }

        @Override
        public Object call(final Object... args) {
            final target = ((DataflowProcessor) getDelegate())
            /*
             * Explaining the following code:
             *
             * - 'out' holds the list of all input values which need to be forwarded (bound) to output as many
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

            def out = args[0..-2]
            def itr = indexes.collect { args[it] }
            List<List> cmb = itr.combinations()

            for( int i=0; i<cmb.size(); i++ ) {
                List entries = cmb[i]
                int count = 0
                for( int j=0; j<len; j++ ) {
                    if( j in this.indexes ) {
                        out[j] = entries[count++]
                    }
                }

                target.bindAllOutputValues( out as Object[] )
            }

            return null
        }

        @Override
        public Object call(final Object arguments) {
            throw new UnsupportedOperationException()
        }

        @Override
        public Object call() {
            throw new UnsupportedOperationException()
        }
    }

    /**
     * Adapter closure to call the {@link #invokeTask(java.lang.Object)} method
     */
    static class InvokeTaskAdapter extends Closure {

        private int numOfParams

        private TaskProcessor processor

        InvokeTaskAdapter(TaskProcessor p, int n) {
            super(null, null);
            processor = p
            numOfParams = n
        }

        @Override
        public int getMaximumNumberOfParameters() {
            numOfParams
        }

        @Override
        public Class[] getParameterTypes() {
            def result = new Class[numOfParams]
            for( int i=0; i<result.size(); i++ )
                result[i] = Object
            return result
        }

        @Override
        public Object call(final Object arguments) {
            processor.invokeTask(arguments)
            return null
        }

        @Override
        public Object call(final Object... args) {
            processor.invokeTask(args as List)
            return null
        }

        @Override
        public Object call() {
            throw new UnsupportedOperationException()
        }
    }

    static enum RunType {
        SUBMIT('Submitted process'),
        RETRY('Re-submitted process')

        String message;

        RunType(String str) { message=str };
    }

    static final public String TASK_CONTEXT_PROPERTY_NAME = 'task'

    final private static Pattern ENV_VAR_NAME = ~/[a-zA-Z_]+[a-zA-Z0-9_]*/

    final private static Pattern QUESTION_MARK = ~/(\?+)/

    /**
     * Keeps track of the task instance executed by the current thread
     */
    protected final ThreadLocal<TaskRun> currentTask = new ThreadLocal<>()

    /**
     * Unique task index number (run)
     */
    final protected indexCount = new AtomicInteger()

    /**
     * The current workflow execution session
     */
    protected Session session

    /**
     * The script object which defines this task
     */
    protected BaseScript ownerScript

    /**
     * Gpars thread pool
     */
    protected PGroup group = Dataflow.retrieveCurrentDFPGroup()

    /**
     * The processor descriptive name
     */
    protected String name

    /**
     * The piece of code to be execute provided by the user
     */
    protected TaskBody taskBody

    /**
     * The corresponding {@code DataflowProcessor} which will receive and
     * manage accordingly the task inputs
     *
     * note: it must be declared volatile -- issue #41
     */
    protected volatile DataflowProcessor operator

    /**
     * The underlying executor which will run the task
     */
    protected Executor executor

    /**
     * The corresponding task configuration properties, it holds the inputs/outputs
     * definition as well as other execution meta-declaration
     */
    protected ProcessConfig config

    /**
     * Count the number of time an error occurred
     */
    private volatile int errorCount

    /**
     * Lock the work dir creation process
     */
    private static Lock lockWorkDirCreation = new ReentrantLock()

    /**
     * Set to true the very first time the error is shown.
     *
     * Note: it is declared static because the error must be shown only the
     * very first time  for all processes
     */
    private static final AtomicBoolean errorShown = new AtomicBoolean()

    /**
     * Flag set {@code true} when the processor termination has been invoked
     *
     * See {@code #checkProcessTermination}
     */
    protected volatile boolean completed

    protected boolean allScalarValues

    protected boolean hasEachParams

    /**
     * Whenever the process execution is required to be blocking in order to handle
     * shared object in a thread safe manner
     */
    protected boolean blocking

    /**
     * The state is maintained by using an agent
     */
    protected Agent<StateObj> state

    /**
     * Groovy engine used to evaluate dynamic code
     */
    protected Grengine grengine

    /**
     * Whenever the process is executed only once
     */
    protected boolean singleton

    /**
     * Track the status of input ports. When 1 the port is open (waiting for data),
     * when 0 the port is closed (ie. received the STOP signal)
     */
    protected AtomicIntegerArray openPorts

    /**
     * Process ID number. The first is 1, the second 2 and so on ..
     */
    private final int id

    private static int processCount

    /*
     * Initialise the process ID
     *
     * Note: processes are create in a sequential manner (by the main thread that parse the script)
     * so it does not require a synchronized block
     */
    {
        id = ++processCount
        grengine = session && session.classLoader ? new Grengine(session.classLoader) : new Grengine()
    }

    /* for testing purpose - do not remove */
    protected TaskProcessor() { }

    /**
     * Create and initialize the processor object
     *
     * @param name
     * @param executor
     * @param session
     * @param script
     * @param config
     * @param taskBody
     */
    TaskProcessor( String name, Executor executor, Session session, BaseScript script, ProcessConfig config, TaskBody taskBody ) {
        assert executor
        assert session
        assert script
        assert taskBody

        this.executor = executor
        this.session = session
        this.ownerScript = script
        this.config = config
        this.taskBody = taskBody
        this.name = name

    }

    /**
     * @return The processor unique id
     */
    int getId() { id }
  
    /**
     * @return The {@code TaskConfig} object holding the task configuration properties
     */
    ProcessConfig getConfig() { config }

    /**
     * @return The current {@code Session} instance
     */
    Session getSession() { session }

    /**
     * @return The processor name
     */
    String getName() { name }

    /**
     * @return The {@link Executor} associated to this processor
     */
    Executor getExecutor() { executor }

    /**
     * @return The {@code DataflowOperator} underlying this process
     */
    DataflowProcessor getOperator() { operator }

    /**
     * @return The {@code BaseScript} object which represents pipeline script
     */
    BaseScript getOwnerScript() { ownerScript }

    /**
     *  Define the type of script hold by the {@code #code} property
     */
    protected ScriptType getScriptType() { taskBody.type }

    /**
     * @return The user provided script block
     */
    public TaskBody getTaskBody() { taskBody }

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
    def run() {

        // -- check that the task has a body
        if ( !taskBody )
            throw new IllegalStateException("Missing task body for process `$name`")

        // -- check that input set defines at least two elements
        def invalidInputSet = config.getInputs().find { it instanceof SetInParam && it.inner.size()<2 }
        if( invalidInputSet )
            log.warn "Input `set` must define at least two component -- Check process `$name`"

        // -- check that output set defines at least two elements
        def invalidOutputSet = config.getOutputs().find { it instanceof SetOutParam && it.inner.size()<2 }
        if( invalidOutputSet )
            log.warn "Output `set` must define at least two component -- Check process `$name`"

        /**
         * Verify if this process run only one time
         */
        allScalarValues = config.getInputs().allScalarInputs()
        hasEachParams = config.getInputs().any { it instanceof EachInParam }

        /*
         * Normalize input channels
         */
        config.fakeInput()

        /*
         * Normalize the output
         * - even though the output may be empty, let return the stdout as output by default
         */
        if ( config.getOutputs().size() == 0 ) {
            config.fakeOutput()
        }

        // the state agent
        state = new Agent<>(new StateObj(name))
        state.addListener { StateObj old, StateObj obj ->
            try {
                log.trace "<$name> Process state changed to: $obj -- finished: ${obj.isFinished()}"
                if( !completed && obj.isFinished() ) {
                    terminateProcess()
                    completed = true
                }
            }
            catch( Throwable e ) {
                session.abort(e)
            }
        }

        // register the processor
        // note: register the task *before* creating (and starting the dataflow operator) in order
        // a race condition on the processes barrier - this fix issue #43
        session.processRegister(this)

        // create the underlying dataflow operator
        createOperator()

        session.notifyProcessCreate(this)

        /*
         * When there is a single output channel, return let returns that item
         * otherwise return the list
         */
        def result = config.getOutputs().channels
        return result.size() == 1 ? result[0] : result
    }

    /**
     * Template method which extending classes have to override in order to
     * create the underlying *dataflow* operator associated with this processor
     *
     * See {@code DataflowProcessor}
     */

    protected void createOperator() {
        def opInputs = new ArrayList(config.getInputs().getChannels())

        /*
         * check if there are some iterators declaration
         * the list holds the index in the list of all *inputs* for the {@code each} declaration
         */
        List<Integer> iteratorIndexes = []
        config.getInputs().eachWithIndex { param, index ->
            if( param instanceof EachInParam ) {
                log.trace "Process ${name} > got each param: ${param.name} at index: ${index} -- ${param.dump()}"
                iteratorIndexes << index
            }
        }

        /*
         * When one (or more) {@code each} are declared as input, it is created an extra
         * operator which will receive the inputs from the channel (excepts the values over iterate)
         *
         * The operator will *expand* the received inputs, iterating over the user provided value and
         * forwarding the final values the the second *parallel* processor executing the user specified task
         */
        if( iteratorIndexes ) {
            log.debug "Creating *combiner* operator for each param(s) at index(es): ${iteratorIndexes}"

            // don't care about the last channel, being the control channel it doesn't bring real values
            final size = opInputs.size()-1

            // the iterator operator needs to executed just one time
            // thus add a dataflow queue binding a single value and then a stop signal
            def termination = new DataflowQueue<>()
            termination << Boolean.TRUE
            opInputs[size] = termination

            // the channel forwarding the data from the *iterator* process to the target task
            final linkingChannels = new ArrayList(size)
            size.times { linkingChannels[it] = new DataflowQueue() }

            // the script implementing the iterating process
            final forwarder = new ForwardClosure(size, iteratorIndexes)

            // instantiate the iteration process
            def stopAfterFirstRun = allScalarValues
            def interceptor = new BaseProcessInterceptor(opInputs, stopAfterFirstRun)
            def params = [inputs: opInputs, outputs: linkingChannels, maxForks: 1, listeners: [interceptor]]
            session.allOperators << (operator = new DataflowOperator(group, params, forwarder))
            // fix issue #41
            operator.start()

            // set as next inputs the result channels of the iteration process
            // adding the 'control' channel removed previously
            opInputs = new ArrayList(size+1)
            opInputs.addAll( linkingChannels )
            opInputs.add( config.getInputs().getChannels().last() )
        }


        /*
         * define the max forks attribute:
         * - by default the process execution is parallel using the poolSize value
         * - otherwise use the value defined by the user via 'taskConfig'
         */
        def maxForks = session.poolSize
        if( config.maxForks ) {
            maxForks = config.maxForks
            blocking = true
        }
        log.debug "Creating operator > $name -- maxForks: $maxForks"

        /*
         * finally create the operator
         */
        // note: do not specify the output channels in the operator declaration
        // this allows us to manage them independently from the operator life-cycle
        this.singleton = allScalarValues && !hasEachParams
        this.openPorts = createPortsArray(opInputs.size())
        config.getOutputs().setSingleton(singleton)
        def interceptor = new TaskProcessorInterceptor(opInputs, singleton)
        def params = [inputs: opInputs, maxForks: maxForks, listeners: [interceptor] ]
        def invoke = new InvokeTaskAdapter(this, opInputs.size())
        session.allOperators << (operator = new DataflowOperator(group, params, invoke))

        // notify the creation of a new vertex the execution DAG
        NodeMarker.addProcessNode(this, config.getInputs(), config.getOutputs())

        // fix issue #41
        operator.start()

    }

    private AtomicIntegerArray createPortsArray(int size) {
        def result = new AtomicIntegerArray(size)
        for( int i=0; i<size; i++ )
            result.set(i, 1)
        return result
    }


    /**
     * The processor execution body
     *
     * @param processor
     * @param values
     */
    final protected void invokeTask( def args ) {

        // create and initialize the task instance to be executed
        final List values = args instanceof List ? args : [args]
        log.trace "Invoking task > $name"

        // -- create the task run instance
        final task = createTaskRun()
        // -- set the task instance as the current in this thread
        currentTask.set(task)

        // -- validate input lengths
        validateInputSets(values)

        // -- map the inputs to a map and use to delegate closure values interpolation
        final secondPass = [:]
        int count = makeTaskContextStage1(task, secondPass, values)
        makeTaskContextStage2(task, secondPass, count)

        // verify that `when` guard, when specified, is satisfied
        if( !checkWhenGuard(task) )
            return

        // -- resolve the task command script
        task.resolve(taskBody)

        // -- verify if exists a stored result for this case,
        //    if true skip the execution and return the stored data
        if( checkStoredOutput(task) )
            return

        def hash = createTaskHashKey(task)
        checkCachedOrLaunchTask(task, hash, resumable)
    }

    @Memoized
    private List<SetInParam> getDeclaredInputSet() {
        getConfig().getInputs().ofType(SetInParam)
    }

    protected void validateInputSets( List values ) {

        def declaredSets = getDeclaredInputSet()
        for( int i=0; i<declaredSets.size(); i++ ) {
            final param = declaredSets[i]
            final entry = values[param.index]
            final expected = param.inner.size()
            final actual = entry instanceof Collection ? entry.size() : (entry instanceof Map ? entry.size() : 1)

            if( actual != expected ) {
                log.warn1("Input tuple does not match input set cardinality declared by process `$name` -- offending value: $entry", firstOnly: true, cacheKey: this)
            }
        }
    }


    /**
     * @return A string 'she-bang' formatted to the added on top script to be executed.
     * The interpreter to be used define bu the *taskConfig* property {@code shell}
     */
    static shebangLine(shell) {
        assert shell, "Missing 'shell' property in process configuration"

        String result = shell instanceof List ? shell.join(' ') : shell
        if( result.startsWith('/') ) {
            result = '#!' + result
        }
        else {
            result= '#!/usr/bin/env ' + result
        }

        return result

    }

    /**
     * Remove extra leading, trailing whitespace and newlines chars,
     * also if the script does not start with a {@code shebang} line,
     * add the default by using the current {@code #shell} attribute
     */
    static String normalizeScript(String script, shell) {
        assert script != null

        def result = new StringBuilder()
        result << script.stripIndent().trim()
        result << '\n'

        if( result[0] != '#' || result[1] != '!') {
            result.insert(0, shebangLine(shell) +'\n')
        }

        return result.toString()
    }

    /**
     * Given the task script extract the top *she-bang* interpreter declaration removing the {@code #!} characters.
     * @param script The script to be executed
     *
     * @return The interpreter as defined in the she-bang declaration, for example {@code /usr/bin/env perl}
     */
    static String fetchInterpreter( String script ) {
        assert script != null

        if( script[0] == '#' && script[1] == '!') {
            return script.readLines()[0].substring(2)
        }

        return null
    }

    /**
     * Create a new {@code TaskRun} instance, initializing the following properties :
     * <li>{@code TaskRun#id}
     * <li>{@code TaskRun#status}
     * <li>{@code TaskRun#index}
     * <li>{@code TaskRun#name}
     * <li>{@code TaskRun#process}
     *
     * @return The new newly created {@code TaskRun}
     */

    final protected TaskRun createTaskRun() {
        log.trace "Creating a new task > $name"

        def index = indexCount.incrementAndGet()
        def task = new TaskRun(
                id: TaskId.next(),
                index: index,
                processor: this,
                type: scriptType,
                config: config.createTaskConfig(),
                context: new TaskContext(this)
        )

        // setup config
        task.config.index = task.index
        task.config.process = task.processor.name
        task.config.executor = task.processor.executor.name

        /*
         * initialize the inputs/outputs for this task instance
         */
        config.getInputs().each { InParam param ->
            if( param instanceof SetInParam )
                param.inner.each { task.setInput(it)  }
            else if( param instanceof EachInParam )
                task.setInput(param.inner)
            else
                task.setInput(param)
        }

        config.getOutputs().each { OutParam param ->
            if( param instanceof SetOutParam ) {
                param.inner.each { task.setOutput(it) }
            }
            else
                task.setOutput(param)
        }

        return task
    }

    /**
     * Try to check if exists a previously executed process result in the a cached folder. If it exists
     * use the that result and skip the process execution, otherwise the task is sumitted for execution.
     *
     * @param task
     *      The {@code TaskRun} instance to be executed
     * @param hash
     *      The unique {@code HashCode} for the given task inputs
     * @param script
     *      The script to be run (only when it's a merge task)
     * @return
     *      {@code false} when a cached result has been found and the execution has skipped,
     *      or {@code true} if the task has been submitted for execution
     *
     */

    final protected boolean checkCachedOrLaunchTask( TaskRun task, HashCode hash, boolean shouldTryCache ) {

        int tries = task.failCount +1
        while( true ) {
            hash = CacheHelper.defaultHasher().newHasher().putBytes(hash.asBytes()).putInt(tries).hash()

            boolean exists=false
            final folder = FileHelper.getWorkFolder(session.workDir, hash)
            lockWorkDirCreation.lock()
            try {
                exists = folder.exists()
                if( !exists && !folder.mkdirs() )
                    throw new IOException("Unable to create folder=$folder -- check file system permission")
            }
            finally {
                lockWorkDirCreation.unlock()
            }

            log.trace "[${task.name}] Cacheable folder=$folder -- exists=$exists; try=$tries; shouldTryCache=$shouldTryCache"
            def cached = shouldTryCache && exists && checkCachedOutput(task.clone(), folder, hash)
            if( cached )
                return false

            if( exists ) {
                tries++
                continue
            }

            // submit task for execution
            submitTask( task, hash, folder )
            return true
        }

    }

    /**
     * Check if exists a *storeDir* for the specified task. When if exists
     * and contains the expected result files, the process execution is skipped.
     *
     * @param task The task for which check the stored output
     * @return {@code true} when the folder exists and it contains the expected outputs,
     *      {@code false} otherwise
     */
    final boolean checkStoredOutput( TaskRun task ) {
        if( !task.config.storeDir ) {
            log.trace "[$task.name] Store dir not set -- return false"
            return false
        }

        // -- when store path is set, only output params of type 'file' can be specified
        final ctx = task.context
        def invalid = task.getOutputs().keySet().any {
            if( it instanceof ValueOutParam ) {
                return !ctx.containsKey(it.name)
            }
            if( it instanceof FileOutParam ) {
                return false
            }
            return true
        }
        if( invalid ) {
            log.warn "[$task.name] StoreDir can only be used when using 'file' outputs"
            return false
        }

        if( !task.config.getStoreDir().exists() ) {
            log.trace "[$task.name] Store dir does not exists > ${task.config.storeDir} -- return false"
            // no folder -> no cached result
            return false
        }


        try {
            final exit = task.config.getValidExitStatus()[0]
            // -- expose task exit status to make accessible as output value
            task.config.exitStatus = exit
            // -- check if all output resources are available
            collectOutputs(task)
            log.info "[skipping] Stored process > ${task.name}"

            // set the exit code in to the task object
            task.exitStatus = exit
            task.cached = true

            // -- now bind the results
            finalizeTask0(task)
            return true
        }
        catch( MissingFileException | MissingValueException e ) {
            log.trace "[$task.name] Missed store > ${e.getMessage()} -- folder: ${task.config.storeDir}"
            task.exitStatus = Integer.MAX_VALUE
            task.workDir = null
            return false
        }
    }

    /**
     * Check whenever the outputs for the specified task already exist
     *
     * @param task The task instance
     * @param folder The folder where the outputs are stored (eventually)
     * @return {@code true} when all outputs are available, {@code false} otherwise
     */
    final boolean checkCachedOutput(TaskRun task, Path folder, HashCode hash) {

        // check if exists the task exit code file
        def exitCode = null
        def exitFile = folder.resolve(TaskRun.CMD_EXIT)
        if( task.type == ScriptType.SCRIPTLET ) {
            def str
            try {
                str = exitFile.text?.trim()
            }
            catch( IOException e ) {
                log.trace "[$task.name] Exit file can't be read > $exitFile -- return false -- Cause: ${e.message}"
                return false
            }

            exitCode = str.isInteger() ? str.toInteger() : null
            if( !task.isSuccess(exitCode) ) {
                log.trace "[$task.name] Exit code is not valid > $str -- return false"
                return false
            }
        }

        /*
         * load the task record in the cache DB
         */

        /*
         * verify cached context map
         */
        TaskEntry entry
        try {
            entry = session.cache.getTaskEntry(hash, task.processor)
            if( !entry ) {
                log.trace "[$task.name] Missing cache entry -- return false"
                return false
            }

            if( task.hasCacheableValues() && !entry.context ) {
                log.trace "[$task.name] Missing cache context -- return false"
                return false
            }

        }
        catch( Throwable e ) {
            log.warn1("[$task.name] Unable to resume cached task -- See log file for details", causedBy: e)
            return false
        }

        /*
         * verify stdout file
         */
        final stdoutFile = folder.resolve( TaskRun.CMD_OUTFILE )

        if( entry.context != null ) {
            task.context = entry.context
            task.config.context = entry.context
            task.code?.delegate = entry.context
        }

        try {
            // -- expose task exit status to make accessible as output value
            task.config.exitStatus = exitCode
            // -- check if all output resources are available
            collectOutputs(task, folder, stdoutFile, task.context)

            // set the exit code in to the task object
            task.cached = true
            task.hash = hash
            task.workDir = folder
            task.stdout = stdoutFile
            if( exitCode != null ) {
                task.exitStatus = exitCode
            }

            log.info "[${task.hashLog}] Cached process > ${task.name}"
            // -- notify cached event
            if( entry )
                session.notifyTaskCached(new CachedTaskHandler(task,entry.trace))

            // -- now bind the results
            finalizeTask0(task)
            return true
        }
        catch( MissingFileException | MissingValueException e ) {
            log.trace "[$task.name] Missed cache > ${e.getMessage()} -- folder: $folder"
            task.exitStatus = Integer.MAX_VALUE
            task.workDir = null
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
        log.trace "Handling error: $error -- task: $task"
        def fault = resumeOrDie(task, error)

        if (fault instanceof TaskFault) {
            session.fault(fault)
            // when a `TaskFault` is returned a `TERMINATE` is implicit, thus return `true`
            return true
        }

        return fault == TERMINATE || fault == FINISH
    }

    /**
     * @param task The {@code TaskRun} instance that raised an error
     * @param error The error object
     * @return
     *      Either a value of value of {@link ErrorStrategy} representing the error strategy chosen
     *      or an instance of {@TaskFault} representing the cause of the error (that implicitly means
     *      a {@link ErrorStrategy#TERMINATE})
     */
    @PackageScope
    final synchronized resumeOrDie( TaskRun task, Throwable error ) {
        log.trace "Handling unexpected condition for\n  task: $task\n  error [${error?.class?.name}]: ${error?.getMessage()?:error}"

        ErrorStrategy errorStrategy = TERMINATE
        final message = []
        try {
            // -- do not recoverable error, just re-throw it
            if( error instanceof Error ) throw error

            // -- retry without increasing the error counts
            if( task && error.cause instanceof CloudSpotTerminationException ) {
                log.info "[$task.hashLog] NOTE: ${error.message} -- Cause: ${error.cause.message} -- Execution is retried"
                task.failCount+=1
                final taskCopy = task.makeCopy()
                taskCopy.runType = RunType.RETRY
                session.getExecService().submit { checkCachedOrLaunchTask( taskCopy, taskCopy.hash, false ) }
                task.failed = true
                task.errorAction = RETRY
                return RETRY
            }

            final int taskErrCount = task ? ++task.failCount : 0
            final int procErrCount = ++errorCount

            // -- when is a task level error and the user has chosen to ignore error,
            //    just report and error message and DO NOT stop the execution
            if( task && error instanceof ProcessException ) {
                // expose current task exist status
                task.config.exitStatus = task.exitStatus
                task.config.errorCount = procErrCount
                task.config.retryCount = taskErrCount

                errorStrategy = checkErrorStrategy(task, error, taskErrCount, procErrCount)
                if( errorStrategy.soft ) {
                    def msg = "[$task.hashLog] NOTE: $error.message"
                    if( errorStrategy == IGNORE ) msg += " -- Error is ignored"
                    else if( errorStrategy == RETRY ) msg += " -- Execution is retried ($taskErrCount)"
                    log.info msg
                    task.failed = true
                    task.errorAction = errorStrategy
                    return errorStrategy
                }
            }

            // -- mark the task as failed
            if( task ) {
                task.failed = true
                task.errorAction = errorStrategy
            }

            // -- make sure the error is showed only the very first time across all processes
            if( errorShown.getAndSet(true) || session.aborted ) {
                return errorStrategy
            }

            def dumpStackTrace = log.isTraceEnabled()
            message << "Error executing process > '${task?.name ?: name}'"
            switch( error ) {
                case ProcessException:
                    formatTaskError( message, error, task )
                    break

                case FailedGuardException:
                    formatGuardError( message, error as FailedGuardException, task )
                    break;

                default:
                    message << formatErrorCause(error)
                    dumpStackTrace = true
            }
            if( dumpStackTrace )
                log.error(message.join('\n'), error)
            else
                log.error(message.join('\n'))
        }
        catch( Throwable e ) {
            // no recoverable error
            log.error("Execution aborted due to an unexpected error", e )
        }

        return new TaskFault(error: error, task: task, report: message.join('\n'))
    }

    protected ErrorStrategy checkErrorStrategy( TaskRun task, ProcessException error, final int taskErrCount, final int procErrCount ) {

        final action = task.config.getErrorStrategy()

        // retry is not allowed when the script cannot be compiled or similar errors
        if( error instanceof ProcessUnrecoverableException ) {
            return !action.soft ? action : TERMINATE
        }

        // IGNORE strategy -- just continue
        if( action == IGNORE ) {
            return IGNORE
        }

        // RETRY strategy -- check that process do not exceed 'maxError' and the task do not exceed 'maxRetries'
        if( action == RETRY ) {
            final int maxErrors = task.config.getMaxErrors()
            final int maxRetries = task.config.getMaxRetries()

            if( (procErrCount < maxErrors || maxErrors == -1) && taskErrCount <= maxRetries ) {
                final taskCopy = task.makeCopy()
                session.getExecService().submit({
                    try {
                        taskCopy.config.attempt = taskErrCount+1
                        taskCopy.runType = RunType.RETRY
                        taskCopy.resolve(taskBody)
                        checkCachedOrLaunchTask( taskCopy, taskCopy.hash, false )
                    }
                    catch( Throwable e ) {
                        log.error("Unable to re-submit task `${taskCopy.name}`", e)
                        session.abort(e)
                    }
                } as Runnable)
                return RETRY
            }

            return TERMINATE
        }

        return action
    }

    final protected List<String> formatGuardError( List<String> message, FailedGuardException error, TaskRun task ) {
        // compose a readable error message
        message << formatErrorCause(error)

        if( error.source )  {
            message << "\nWhen block:"
            error.source.stripIndent().eachLine {
                message << "  $it"
            }
        }

        if( task?.workDir )
            message << "\nWork dir:\n  ${task.workDirStr}"

        return message
    }

    final protected List<String> formatTaskError( List<String> message, Throwable error, TaskRun task ) {

        // compose a readable error message
        message << formatErrorCause( error )

        /*
         * task executing scriptlets
         */
        if( task?.script ) {
            // - print the executed command
            message << "Command executed${task.template ? " [$task.template]": ''}:\n"
            task.script?.stripIndent()?.trim()?.eachLine {
                message << "  ${it}"
            }

            // - the exit status
            message << "\nCommand exit status:\n  ${task.exitStatus != Integer.MAX_VALUE ? task.exitStatus : '-'}"

            // - the tail of the process stdout
            message << "\nCommand output:"
            final max = 50
            def lines = task.dumpStdout(max)
            if( lines.size() == 0 ) {
                message << "  (empty)"
            }
            lines.each {
                message << "  ${task.workDir ? it.replace(task.workDir.toString()+'/','') : it }"
            }

            // - the tail of the process stderr
            lines = task.dumpStderr(max)
            if( lines ) {
                message << "\nCommand error:"
                lines.each {
                    message << "  ${task.workDir ? it.replace(task.workDir.toString()+'/','') : it }"
                }
            }
            // - this is likely a task wrapper issue
            else if( task.exitStatus != 0 ) {
                lines = task.dumpLogFile(max)
                if( lines ) {
                    message << "\nCommand wrapper:"
                    lines.each {
                        message << "  ${task.workDir ? it.replace(task.workDir.toString()+'/','') : it }"
                    }
                }
            }

        }
        else {
            if( task?.source )  {
                message << "Source block:"
                task.source.stripIndent().eachLine {
                    message << "  $it"
                }
            }

        }

        if( task?.workDir )
            message << "\nWork dir:\n  ${task.workDirStr}"

        message << "\nTip: ${getRndTip()}"

        return message
    }

    static List tips = [
            'when you have fixed the problem you can continue the execution appending to the nextflow command line the option `-resume`',
            "you can try to figure out what's wrong by changing to the process work dir and showing the script file named `${TaskRun.CMD_SCRIPT}`",
            "view the complete command output by changing to the process work dir and entering the command `cat ${TaskRun.CMD_OUTFILE}`",
            "you can replicate the issue by changing to the process work dir and entering the command `bash ${TaskRun.CMD_RUN}`"
    ]

    static Random RND = Random.newInstance()

    /**
     * Display a random tip at the bottom of the error report
     *
     * @return The tip string to display
     */
    protected String getRndTip() {
        tips[ RND.nextInt( tips.size() ) ]
    }


    /**
     * Send a poison pill over all the outputs channel
     */
    final protected synchronized void sendPoisonPill() {
        log.trace "<$name> Sending a poison pill(s)"

        config.getOutputs().getChannels().each { channel ->

            if( channel instanceof DataflowQueue ) {
                channel.bind( PoisonPill.instance )
            }
            else if( channel instanceof DataflowStreamWriteAdapter ) {
                channel.bind( PoisonPill.instance )
            }
            else if( channel instanceof DataflowExpression && !channel.isBound()) {
                channel.bind( PoisonPill.instance )
            }
        }
    }


    private String formatErrorCause( Throwable error ) {

        def result = new StringBuilder()
        result << '\nCaused by:\n'

        def message
        if( error instanceof ProcessStageException || error instanceof MissingFileException || !error.cause )
            message = error.getMessage()
        else
            message = error.cause.getMessage()

        if( !message )
            message = error.toString()

        result
            .append('  ')
            .append(message)
            .append('\n')
            .toString()
    }

    /**
     * Publish output files to a specified target folder
     *
     * @param task The task whose outputs need to be published
     * @param overwrite When {@code true} any existing file will be overwritten, otherwise the publishing is ignored
     */
    @CompileStatic
    protected void publishOutputs( TaskRun task ) {
        def publish = task.config.getPublishDir()
        if( !publish ) {
            return
        }

        if( publish.overwrite == null ) {
            publish.overwrite = !task.cached
        }

        List<Path> files = []
        task.getOutputsByType(FileOutParam).each { param, value ->
            if( value instanceof Path ) {
                files << ((Path)value)
            }
            else if( value instanceof Collection<Path> ) {
                value.each { Path it -> files << it }
            }
            else if( value != null ) {
                throw new IllegalArgumentException("Unknown output file object [${value.class.name}]: ${value}")
            }
        }


        publish.apply(files, task)
    }


    /**
     * Bind the expected output files to the corresponding output channels
     * @param processor
     */
    synchronized protected void bindOutputs( TaskRun task ) {

        // -- creates the map of all tuple values to bind
        Map<Short,List> tuples = [:]
        for( OutParam param : config.getOutputs() ) {
            tuples.put(param.index, [])
        }

        // -- collects the values to bind
        for( OutParam param: task.outputs.keySet() ){
            def value = task.outputs.get(param)

            switch( param ) {
            case StdOutParam:
                log.trace "Process $name > normalize stdout param: $param"
                value = value instanceof Path ? value.text : value?.toString()

            case OptionalParam:
                if( !value && param instanceof OptionalParam && param.optional ) {
                    final holder = [] as MissingParam; holder.missing = param
                    tuples[param.index] = holder
                    break
                }

            case ValueOutParam:
                log.trace "Process $name > collecting out param: ${param} = $value"
                tuples[param.index].add(value)
                break

            default:
                throw new IllegalArgumentException("Illegal output parameter type: $param")
            }
        }

        // -- bind out the collected values
        for( OutParam param : config.getOutputs() ) {
            def list = tuples[param.index]
            if( list == null )
                throw new IllegalStateException()

            if( list instanceof MissingParam ) {
                log.debug "Process $name > Skipping output binding because one or more optional files are missing: $list.missing"
                continue
            }

            if( param.mode == BasicMode.standard ) {
                log.trace "Process $name > Binding out param: ${param} = ${list}"
                bindOutParam(param, list)
            }

            else if( param.mode == BasicMode.flatten ) {
                log.trace "Process $name > Flatting out param: ${param} = ${list}"
                CollectionHelper.flatten( list ) {
                    bindOutParam( param, it )
                }
            }

            else if( param.mode == SetOutParam.CombineMode.combine ) {
                log.trace "Process $name > Combining out param: ${param} = ${list}"
                def combs = list.combinations()
                combs.each { bindOutParam(param, it) }
            }

            else
                throw new IllegalStateException("Unknown bind output parameter type: ${param}")
        }

        // -- finally prints out the task output when 'echo' is true
        if( task.config.echo ) {
            task.echoStdout()
        }
    }

    protected void bindOutParam( OutParam param, List values ) {
        log.trace "<$name> Binding param $param with $values"
        def x = values.size() == 1 ? values[0] : values
        param.getOutChannels().each { it.bind(x) }
    }

    protected void collectOutputs( TaskRun task ) {
        collectOutputs( task, task.getTargetDir(), task.@stdout, task.context )
    }

    /**
     * Once the task has completed this method is invoked to collected all the task results
     *
     * @param task
     */
    final protected void collectOutputs( TaskRun task, Path workDir, def stdout, Map context ) {
        log.trace "<$name> collecting output: ${task.outputs}"

        task.outputs.keySet().each { OutParam param ->

            switch( param ) {
                case StdOutParam:
                    collectStdOut(task, (StdOutParam)param, stdout)
                    break

                case FileOutParam:
                    collectOutFiles(task, (FileOutParam)param, workDir, context)
                    break

                case ValueOutParam:
                    collectOutValues(task, (ValueOutParam)param, context)
                    break

                default:
                    throw new IllegalArgumentException("Illegal output parameter: ${param.class.simpleName}")

            }
        }

        // mark ready for output binding
        task.canBind = true
    }

    /**
     * Collects the process 'std output'
     *
     * @param task The executed process instance
     * @param param The declared {@link StdOutParam} object
     * @param stdout The object holding the task produced std out object
     */
    protected void collectStdOut( TaskRun task, StdOutParam param, def stdout ) {

        if( stdout == null && task.type == ScriptType.SCRIPTLET ) {
            throw new IllegalArgumentException("Missing 'stdout' for process > ${task.name}")
        }

        if( stdout instanceof Path && !stdout.exists() ) {
            throw new MissingFileException("Missing 'stdout' file: ${stdout} for process > ${task.name}")
        }

        task.setOutput(param, stdout)
    }


    protected void collectOutFiles( TaskRun task, FileOutParam param, Path workDir, Map context ) {

        final List<Path> allFiles = []
        // type file parameter can contain a multiple files pattern separating them with a special character
        def entries = param.getFilePatterns(context, task.workDir)

        // for each of them collect the produced files
        entries.each { String filePattern ->
            List<Path> result = null

            def splitter = param.glob ? FilePatternSplitter.glob().parse(filePattern) : null
            if( splitter?.isPattern() ) {
                result = fetchResultFiles(param, filePattern, workDir)
                // filter the inputs
                if( !param.includeInputs ) {
                    result = filterByRemovingStagedInputs(task, result)
                    log.trace "Process ${task.name} > after removing staged inputs: ${result}"
                }
            }
            else {
                def path = param.glob ? splitter.strip(filePattern) : filePattern
                def file = workDir.resolve(path)
                def exists = param.followLinks ? file.exists() : file.exists(LinkOption.NOFOLLOW_LINKS)
                if( exists )
                    result = [file]
                else
                    log.debug "Process `${task.name}` is unable to find [${file.class.simpleName}]: `$file` (pattern: `$filePattern`)"
            }

            if( result )
                allFiles.addAll(result)

            else if( !param.optional )
                throw new MissingFileException("Missing output file(s) `$filePattern` expected by process `${task.name}`")
        }

        task.setOutput( param, allFiles.size()==1 ? allFiles[0] : allFiles )

    }


    protected void collectOutValues( TaskRun task, ValueOutParam param, Map ctx ) {

        try {
            // fetch the output value
            final val = param.resolve(ctx)
            // set into the output set
            task.setOutput(param,val)
            // trace the result
            log.trace "Collecting param: ${param.name}; value: ${val}"
        }
        catch( MissingPropertyException e ) {
            throw new MissingValueException("Missing value declared as output parameter: ${e.property}")
        }

    }

    /**
     * Collect the file(s) with the name specified, produced by the execution
     *
     * @param workDir The job working path
     * @param namePattern The file name, it may include file name wildcards
     * @return The list of files matching the specified name in lexicographical order
     * @throws MissingFileException when no matching file is found
     */
    @PackageScope
    List<Path> fetchResultFiles( FileOutParam param, String namePattern, Path workDir ) {
        assert namePattern
        assert workDir

        List files = []
        def opts = visitOptions(param, namePattern)
        // scan to find the file with that name
        try {
            FileHelper.visitFiles(opts, workDir, namePattern) { Path it -> files.add(it) }
        }
        catch( NoSuchFileException e ) {
            throw new MissingFileException("Cannot access folder: '$workDir'", e)
        }

        return files.sort()
    }

    /**
     * Given a {@link FileOutParam} object create the option map for the
     * {@link FileHelper#visitFiles(java.util.Map, java.nio.file.Path, java.lang.String, groovy.lang.Closure)} method
     *
     * @param param A task {@link FileOutParam}
     * @param namePattern A file glob pattern
     * @return A {@link Map} object holding the traverse options for the {@link FileHelper#visitFiles(java.util.Map, java.nio.file.Path, java.lang.String, groovy.lang.Closure)} method
     */
    @PackageScope
    Map visitOptions( FileOutParam param, String namePattern ) {
        final opts = [:]
        opts.relative = false
        opts.hidden = param.hidden ?: namePattern.startsWith('.')
        opts.followLinks = param.followLinks
        opts.maxDepth = param.maxDepth
        opts.type = param.type ? param.type : ( namePattern.contains('**') ? 'file' : 'any' )
        return opts
    }

    /**
     * Given a list of {@code Path} removes all the hidden file i.e. the ones which names starts with a dot char
     * @param files A list of {@code Path}
     * @return The result list not containing hidden file entries
     */
    @PackageScope
    List<Path> filterByRemovingHiddenFiles( List<Path> files ) {
        files.findAll { !it.getName().startsWith('.') }
    }

    /**
     * Given a list of {@code Path} removes all the entries which name match the name of
     * file used as input for the specified {@code TaskRun}
     *
     * See TaskRun#getStagedInputs
     *
     * @param task
     * @param files
     * @return
     */
    @PackageScope
    List<Path> filterByRemovingStagedInputs( TaskRun task, List<Path> files ) {

        // get the list of input files
        def List<String> allStaged = task.getStagedInputs()
        files.findAll { !allStaged.contains(it.getName()) }

    }


    /**
     * @return The map holding the shell environment variables for the task to be executed
     */
    @Memoized
    def Map<String,String> getProcessEnvironment() {

        def result = [:]

        // add the taskConfig environment entries
        if( session.config.env instanceof Map ) {
            session.config.env.each { name, value ->
                result.put( name, value?.toString() )
            }
        }
        else {
            log.debug "Invalid 'session.config.env' object: ${session.config.env?.class?.name}"
        }


        // pre-pend the 'bin' folder to the task environment
        if( session.binDir ) {
            if( result.containsKey('PATH') ) {
                // note: do not escape potential blanks in the bin path because the PATH
                // variable is enclosed in `"` when in rendered in the launcher script -- see #630
                result['PATH'] =  "${session.binDir}:${result['PATH']}".toString()
            }
            else {
                result['PATH'] = "${session.binDir}:\$PATH".toString()
            }
        }

        return Collections.unmodifiableMap(result)
    }

    /**
     * An input file parameter can be provided with any value other than a file.
     * This function normalize a generic value to a {@code Path} create a temporary file
     * in the for it.
     *
     * @param input The input value
     * @param altName The name to be used when a temporary file is created.
     * @return The {@code Path} that will be staged in the task working folder
     */
    protected FileHolder normalizeInputToFile( Object input, String altName ) {

        /*
         * when it is a local file, just return a reference holder to it
         */
        if( input instanceof Path ) {
            return new FileHolder(input)
        }

        /*
         * default case, convert the input object to a string and save
         * to a local file
         */
        def source = input?.toString() ?: ''
        def result = Nextflow.tempFile(altName)
        result.text = source
        return new FileHolder(source, result)
    }


    protected List<FileHolder> normalizeInputToFiles( Object obj, int count ) {

        Collection allItems = obj instanceof Collection ? obj : [obj]
        def len = allItems.size()

        // use a bag so that cache hash key is not affected by file entries order
        def files = new ArrayBag<FileHolder>(len)
        for( def item : allItems ) {
            files << normalizeInputToFile(item, "input.${++count}")
        }

        return files
    }

    protected singleItemOrList( List<FileHolder> items ) {
        assert items != null

        if( items.size() == 1 ) {
            return Paths.get(items[0].stageName)
        }

        def result = new ArrayList(items.size())
        for( int i=0; i<items.size(); i++ ) {
            result.add( Paths.get(items[i].stageName) )
        }
        return new BlankSeparatedList(result)
    }


    /**
     * An input file name may contain wildcards characters which have to be handled coherently
     * given the number of files specified.
     *
     * @param name A file name with may contain a wildcard character star {@code *} or question mark {@code ?}.
     *  Only one occurrence can be specified for star or question mark wildcards.
     *
     * @param value Any value that have to be managed as an input files. Values other than {@code Path} are converted
     * to a string value, using the {@code #toString} method and saved in the local file-system. Value of type {@code Collection}
     * are expanded to multiple values accordingly.
     *
     * @return
     */
    @CompileStatic
    protected List<FileHolder> expandWildcards( String name, List<FileHolder> files ) {
        assert files != null

        // use an unordered so that cache hash key is not affected by file entries order
        final result = new ArrayBag(files.size())
        if( files.size()==0 ) { return result }

        if( !name || name == '*' ) {
            result.addAll(files)
            return result
        }

        if( !name.contains('*') && !name.contains('?') && files.size()>1 ) {
            /*
             * The name do not contain any wildcards *BUT* when multiple files are provide
             * it is managed like having a 'star' at the end of the file name
             */
            name += '*'
        }

        for( int i=0; i<files.size(); i++ ) {
            def holder = files[i]
            def newName = expandWildcards0(name, holder.stageName, i+1, files.size())
            result << holder.withName( newName )
        }

        return result
    }

    @CompileStatic
    protected String replaceQuestionMarkWildcards(String name, int index) {
        def result = new StringBuffer()

        Matcher m = QUESTION_MARK.matcher(name)
        while( m.find() ) {
            def match = m.group(1)
            def repString = String.valueOf(index).padLeft(match.size(), '0')
            m.appendReplacement(result, repString)
        }
        m.appendTail(result)
        result.toString()
    }

    @CompileStatic
    protected String replaceStarWildcards(String name, int index, boolean strip=false) {
        name.replaceAll(/\*/, strip ? '' : String.valueOf(index))
    }

    @CompileStatic
    protected String expandWildcards0( String path, String stageName, int index, int size ) {

        String name
        String parent
        int p = path.lastIndexOf('/')
        if( p == -1 ) {
            parent = null
            name = path
        }
        else {
            parent = path.substring(0,p)
            name = path.substring(p+1)
        }

        if( name == '*' || !name ) {
            name = stageName
        }
        else {
            final stripWildcard = size<=1 // <-- string the start wildcard instead of expanding to an index number when the collection contain only one file
            name = replaceStarWildcards(name, index, stripWildcard)
            name = replaceQuestionMarkWildcards(name, index)
        }

        if( parent ) {
            parent = replaceStarWildcards(parent, index)
            parent = replaceQuestionMarkWildcards(parent, index)
            return "$parent/$name"
        }
        else {
            return name
        }

    }

    /**
     * Given a map holding variables key-value pairs, create a script fragment
     * exporting the required environment variables
     */
    static String bashEnvironmentScript( Map<String,String> environment, boolean escape=false ) {
        if( !environment )
            return null

        final List script = []
        environment.each { name, value ->
            if( !ENV_VAR_NAME.matcher(name).matches() )
                log.trace "Illegal environment variable name: '${name}' -- This variable definition is ignored"
            else if( !value ) {
                log.warn "Environment variable `$name` evaluates to an empty value"
                script << "export $name=''"
            }
            else
                script << "export $name=\"${escape ? value.replace('$','\\$') : value}\""
        }
        script << ''

        return script.join('\n')
    }

    final protected int makeTaskContextStage1( TaskRun task, Map secondPass, List values ) {

        final contextMap = task.context
        int count = 0

        task.inputs.keySet().each { InParam param ->

            // add the value to the task instance
            def val = param.decodeInputs(values)

            switch(param) {
                case ValueInParam:
                    contextMap.put( param.name, val )
                    break

                case FileInParam:
                    secondPass[param] = val
                    return // <-- leave it, because we do not want to add this 'val' at this stage

                case StdInParam:
                case EnvInParam:
                    // nothing to do
                    break

                default:
                    throw new IllegalStateException("Unsupported input param type: ${param?.class?.simpleName}")
            }

            // add the value to the task instance context
            task.setInput(param, val)
        }

        return count
    }

    final protected void makeTaskContextStage2( TaskRun task, Map secondPass, int count ) {

        final ctx = task.context
        final allNames = new HashMap<String,Integer>()

        // -- all file parameters are processed in a second pass
        //    so that we can use resolve the variables that eventually are in the file name
        secondPass.each { FileInParam param, val ->

            def fileParam = param as FileInParam
            def normalized = normalizeInputToFiles(val,count)
            def resolved = expandWildcards( fileParam.getFilePattern(ctx), normalized )
            ctx.put( param.name, singleItemOrList(resolved) )
            count += resolved.size()
            for( FileHolder item : resolved ) {
                Integer num = allNames.getOrCreate(item.stageName, 0) +1
                allNames.put(item.stageName,num)
            }

            // add the value to the task instance context
            task.setInput(param, resolved)
        }

        // -- set the delegate map as context ih the task config
        //    so that lazy directives will be resolved against it
        task.config.context = ctx

        // check conflicting file names
        def conflicts = allNames.findAll { name, num -> num>1 }
        if( conflicts ) {
            log.debug("Process $name > collision check staing file names: $allNames")
            def message = "Process `$name` input file name collision -- There are multiple input files for each of the following file names: ${conflicts.keySet().join(', ')}"
            throw new ProcessUnrecoverableException(message)
        }
    }

    final protected void makeTaskContextStage3( TaskRun task, HashCode hash, Path folder ) {

        // set hash-code & working directory
        task.hash = hash
        task.workDir = folder
        task.config.workDir = folder
        task.config.hash = hash.toString()
        task.config.name = task.getName()

    }

    final protected HashCode createTaskHashKey(TaskRun task) {

        List keys = [ session.uniqueId, name, task.source ]

        if( task.containerEnabled )
            keys << task.container

        // add all the input name-value pairs to the key generator
        task.inputs.each {
            keys.add( it.key.name )
            keys.add( it.value )
        }

        // add all variable references in the task script but not declared as input/output
        def vars = getTaskGlobalVars(task)
        if( vars ) {
            log.trace "Task: $name > Adding script vars hash code: ${vars}"
            vars.each { k, v -> keys.add( k ); keys.add( v ) }
        }

        final binEntries = getTaskBinEntries(task.source)
        if( binEntries ) {
            log.trace "Task: $name > Adding scripts on project bin path: ${-> binEntries.join('; ')}"
            keys.addAll(binEntries)
        }

        final modules = task.getConfig().getModule()
        if( modules ) {
            keys.addAll(modules)
        }

        final mode = config.getHashMode()
        final hash = CacheHelper.hasher(keys, mode).hash()
        if( session.dumpHashes ) {
            traceInputsHashes(task, keys, mode, hash)
        }
        return hash
    }

    /**
     * This method scans the task command string looking for invocations of scripts
     * defined in the project bin folder.
     *
     * @param script The task command string
     * @return The list of paths of scripts in the project bin folder referenced in the task command
     */
    @Memoized
    protected List<Path> getTaskBinEntries(String script) {
        def result = []
        def tokenizer = new StringTokenizer(script," \t\n\r\f()[]{};&|<>`")
        while( tokenizer.hasMoreTokens() ) {
            def token = tokenizer.nextToken()
            def path = session.binEntries.get(token)
            if( path )
                result.add(path)
        }
        return result;
    }

    private void traceInputsHashes( TaskRun task, List entries, CacheHelper.HashMode mode, hash ) {

        def buffer = new StringBuilder()
        buffer.append("[${task.name}] cache hash: ${hash}; mode: $mode; entries: \n")
        for( Object item : entries ) {
            buffer.append( "  ${CacheHelper.hasher(item, mode).hash()} [${item?.getClass()?.getName()}] $item \n")
        }

        log.info(buffer.toString())
    }

    protected Map<String,Object> getTaskGlobalVars(TaskRun task) {
        getTaskGlobalVars( task.context.getVariableNames(), ownerScript.binding, task.context.getHolder() )
    }

    /**
     * @param variableNames The collection of variables referenced in the task script
     * @param binding The script global binding
     * @param context The task variable context
     * @return The set of task variables accessed in global script context and not declared as input/output
     */
    final protected Map<String,Object> getTaskGlobalVars(Set<String> variableNames, Binding binding, Map context) {
        final result = new HashMap(variableNames.size())
        final processName = name

        def itr = variableNames.iterator()
        while( itr.hasNext() ) {
            final varName = itr.next()

            final p = varName.indexOf('.')
            final baseName = p !=- 1 ? varName.substring(0,p) : varName

            if( context.containsKey(baseName) ) {
                // when the variable belong to the task local context just ignore it,
                // because it must has been provided as an input parameter
                continue
            }

            if( binding.hasVariable(baseName) ) {
                def value
                try {
                    if( p == -1 ) {
                        value = binding.getVariable(varName)
                    }
                    else {
                        value = grengine.run(varName, binding)
                    }
                }
                catch( MissingPropertyException | NullPointerException e ) {
                    value = null
                    log.trace "Process `${processName}` cannot access global variable `$varName` -- Cause: ${e.message}"
                }

                // value for 'workDir' and 'baseDir' folders are added always as string
                // in order to avoid to invalid the cache key when resuming the execution
                if( varName=='workDir' || varName=='baseDir' )
                    value = value.toString()

                result.put( varName, value )
            }

        }
        return result
    }


    /**
     * Execute the specified task shell script
     *
     * @param script The script string to be execute, e.g. a BASH script
     * @return {@code TaskDef}
     */
    final protected void submitTask( TaskRun task, HashCode hash, Path folder ) {
        log.trace "[${task.name}] actual run folder: ${folder}"

        makeTaskContextStage3(task, hash, folder)

        // add the task to the collection of running tasks
        session.dispatcher.submit(task, blocking)

    }

    protected boolean checkWhenGuard(TaskRun task) {

        try {
            def pass = task.config.getGuard('when')
            if( pass ) {
                return true
            }

            log.trace "Task ${task.name} is not executed because `when` condition is not verified"
            finalizeTask0(task)
            return false
        }
        catch ( FailedGuardException error ) {
            handleException(error, task)
            return false
        }
    }

    /**
     * Finalize the task execution, checking the exit status
     * and binding output values accordingly
     *
     * @param task The {@code TaskRun} instance to finalize
     */
    @PackageScope
    final finalizeTask( TaskRun task ) {
        log.trace "finalizing process > ${task.name} -- $task"

        def fault = null
        try {
            // -- verify task exist status
            if( task.error )
                throw new ProcessFailedException("Process `${task.name}` failed", task.error)

            if( task.type == ScriptType.SCRIPTLET ) {
                if( task.exitStatus == Integer.MAX_VALUE )
                    throw new ProcessFailedException("Process `${task.name}` terminated for an unknown reason -- Likely it has been terminated by the external system")

                if ( !task.isSuccess() )
                    throw new ProcessFailedException("Process `${task.name}` terminated with an error exit status (${task.exitStatus})")
            }

            // -- expose task exit status to make accessible as output value
            task.config.exitStatus = task.exitStatus
            // -- if it's OK collect results and finalize
            collectOutputs(task)
        }
        catch ( Throwable error ) {
            fault = resumeOrDie(task, error)
        }

        // -- finalize the task
        if( fault != ErrorStrategy.RETRY )
            finalizeTask0(task)

        return fault
    }

    /**
     * Whenever the process can be cached
     */
    boolean isCacheable() {
        session.cacheable && config.cacheable
    }

    @PackageScope boolean isResumable() {
        isCacheable() && session.resumeMode
    }

    /**
     * Finalize the task execution, checking the exit status
     * and binding output values accordingly
     *
     * @param task The {@code TaskRun} instance to finalize
     * @param producedFiles The map of files to be bind the outputs
     */
    private void finalizeTask0( TaskRun task ) {
        log.trace "Finalize process > ${task.name}"

        // -- bind output (files)
        if( task.canBind ) {
            bindOutputs(task)
            publishOutputs(task)
        }

        // increment the number of processes executed
        state.update { StateObj it -> it.incCompleted() }
    }

    protected void terminateProcess() {
        log.trace "<${name}> Sending poison pills and terminating process"
        sendPoisonPill()
        session.processDeregister(this)
    }

    /**
     * Dump the current process status listing
     * all input *port* statuses for debugging purpose
     *
     * @return The text description representing the process status
     */
    String dumpTerminationStatus() {

        def result = new StringBuilder()
        def terminated = !DataflowHelper.isProcessorActive(operator)
        result << "[process] $name\n"
        if( terminated )
            return result.toString()

        def statusStr = !completed && !terminated ? 'status=ACTIVE' : ( completed && terminated ? 'status=TERMINATED' : "completed=$completed; terminated=$terminated" )
        result << "  $statusStr\n"
        // add extra info about port statuses
        for( int i=0; i<openPorts.length(); i++ ) {
            def last = i == openPorts.length()-1
            def param = config.getInputs()[i]
            def isValue = param?.inChannel instanceof DataflowExpression
            def type = last ? '(cntrl)' : (isValue ? '(value)' : '(queue)')
            def channel = param && !(param instanceof SetInParam) ? param.getName() : '-'
            def status = type != '(value)' ? (openPorts.get(i) ? 'OPEN' : 'CLOSED') : '-   '
            result << "  port $i: $type ${status}; channel: $channel\n"
        }

        return result.toString()
    }

    /*
     * logger class for the *iterator* processor
     */
    class BaseProcessInterceptor extends DataflowEventAdapter {

        final List<DataflowChannel> inputs

        final boolean stopAfterFirstRun

        final int len

        final DataflowQueue control

        final int first

        BaseProcessInterceptor( List<DataflowChannel> inputs, boolean stop ) {
            this.inputs = new ArrayList<>(inputs)
            this.stopAfterFirstRun = stop
            this.len = inputs.size()
            this.control = (DataflowQueue)inputs.get(len-1)
            this.first = inputs.findIndexOf { it instanceof DataflowQueue }
        }

        @Override
        public Object messageArrived(final DataflowProcessor processor, final DataflowReadChannel<Object> channel, final int index, final Object message) {
            if( len == 1 || stopAfterFirstRun ) {
                // -- kill itself
                control.bind(PoisonPill.instance)
            }
            else if( index == first ) {
                // -- keep things rolling
                control.bind(Boolean.TRUE)
            }

            return message;
        }

    }

    /**
     *  Intercept dataflow process events
     */
    class TaskProcessorInterceptor extends BaseProcessInterceptor {

        TaskProcessorInterceptor(List<DataflowChannel> inputs, boolean stop) {
            super(inputs, stop)
        }

        @Override
        public List<Object> beforeRun(final DataflowProcessor processor, final List<Object> messages) {
            log.trace "<${name}> Before run -- messages: ${messages}"
            // the counter must be incremented here, otherwise it won't be consistent
            state.update { StateObj it -> it.incSubmitted() }
            return messages;
        }


        @Override
        void afterRun(DataflowProcessor processor, List<Object> messages) {
            log.trace "<${name}> After run"
            currentTask.remove()
        }

        @Override
        public Object messageArrived(final DataflowProcessor processor, final DataflowReadChannel<Object> channel, final int index, final Object message) {
            if( log.isTraceEnabled() ) {
                def channelName = config.getInputs()?.names?.get(index)
                def taskName = currentTask.get()?.name ?: name
                log.trace "<${taskName}> Message arrived -- ${channelName} => ${message}"
            }

            super.messageArrived(processor, channel, index, message)
        }

        @Override
        public Object controlMessageArrived(final DataflowProcessor processor, final DataflowReadChannel<Object> channel, final int index, final Object message) {
            if( log.isTraceEnabled() ) {
                def channelName = config.getInputs()?.names?.get(index)
                def taskName = currentTask.get()?.name ?: name
                log.trace "<${taskName}> Control message arrived ${channelName} => ${message}"
            }

            super.controlMessageArrived(processor, channel, index, message)

            if( message == PoisonPill.instance ) {
                log.trace "<${name}> Poison pill arrived; port: $index"
                openPorts.set(index, 0) // mark the port as closed
                state.update { StateObj it -> it.poison() }
            }

            return message
        }

        @Override
        public void afterStop(final DataflowProcessor processor) {
            log.trace "<${name}> After stop"
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
            // return `true` to terminate the dataflow processor
            handleException( error, currentTask.get() )
        }
    }
}
