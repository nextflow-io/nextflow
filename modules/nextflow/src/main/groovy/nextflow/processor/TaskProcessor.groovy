/*
 * Copyright 2013-2026, Seqera Labs
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package nextflow.processor

import java.nio.file.FileSystems
import java.nio.file.Path
import java.util.concurrent.atomic.AtomicInteger
import java.util.concurrent.atomic.AtomicIntegerArray
import java.util.concurrent.atomic.LongAdder
import java.util.function.Consumer
import java.util.regex.Pattern

import ch.artecat.grengine.Grengine
import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import groovyx.gpars.agent.Agent
import groovyx.gpars.dataflow.Dataflow
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.expression.DataflowExpression
import groovyx.gpars.dataflow.operator.DataflowEventAdapter
import groovyx.gpars.dataflow.operator.DataflowOperator
import groovyx.gpars.dataflow.operator.DataflowProcessor
import groovyx.gpars.dataflow.operator.PoisonPill
import groovyx.gpars.dataflow.stream.DataflowStreamWriteAdapter
import groovyx.gpars.group.PGroup
import nextflow.NF
import nextflow.Session
import nextflow.ast.TaskCmdXform
import nextflow.ast.TaskTemplateVarsXform
import nextflow.dag.NodeMarker
import nextflow.exception.ProcessUnrecoverableException
import nextflow.executor.Executor
import nextflow.executor.TaskArrayExecutor
import nextflow.extension.CH
import nextflow.extension.DataflowHelper
import nextflow.script.BaseScript
import nextflow.script.BodyDef
import nextflow.script.ProcessConfig
import nextflow.script.ProcessConfigV1
import nextflow.script.ProcessConfigV2
import nextflow.script.ScriptMeta
import nextflow.script.ScriptType
import nextflow.script.bundle.ResourcesBundle
import nextflow.script.params.DefaultOutParam
import nextflow.script.params.EachInParam
import nextflow.script.params.EnvOutParam
import nextflow.script.params.FileOutParam
import nextflow.script.params.MissingParam
import nextflow.script.params.OptionalParam
import nextflow.script.params.OutParam
import nextflow.script.params.StdOutParam
import nextflow.script.params.TupleInParam
import nextflow.script.params.TupleOutParam
import nextflow.script.params.ValueOutParam
import nextflow.util.Escape
import nextflow.util.TestOnly
import org.codehaus.groovy.control.CompilerConfiguration
import org.codehaus.groovy.control.customizers.ASTTransformationCustomizer
/**
 * Implement nextflow process execution logic
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class TaskProcessor {

    static enum RunType {
        SUBMIT('Submitted process'),
        RETRY('Re-submitted process')

        String message;

        RunType(String str) { message=str };
    }

    static final public String TASK_CONTEXT_PROPERTY_NAME = 'task'

    final private static Pattern ENV_VAR_NAME = ~/[a-zA-Z_]+[a-zA-Z0-9_]*/

    @TestOnly private static volatile TaskProcessor currentProcessor0

    @TestOnly static TaskProcessor currentProcessor() { currentProcessor0 }

    /**
     * Keeps track of the task instance executed by the current thread
     */
    protected final ThreadLocal<TaskRun> currentTask = new ThreadLocal<>()

    /**
     * Unique task index number (run)
     */
    final protected AtomicInteger indexCount = new AtomicInteger()

    /**
     * The current workflow execution session
     */
    protected Session session

    /**
     * The script object which defines this task
     */
    protected BaseScript ownerScript

    /**
     * The processor descriptive name
     */
    protected String name

    /**
     * The piece of code to be execute provided by the user
     */
    protected BodyDef taskBody

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
     * Executes the tasks for this process. It implements the process
     * functionality that does not depend on the dataflow interface, and it
     * reports the completion of each task via {@link #taskCompleted}
     */
    private final TaskRunner runner = new TaskRunner(this, { TaskRun task -> taskCompleted(task) } as Consumer<TaskRun>)

    /**
     * Flag set {@code true} when the processor termination has been invoked
     *
     * See {@code #checkProcessTermination}
     */
    protected volatile boolean completed

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

    private LongAdder forksCount

    private int maxForks

    private static int processCount

    private List<TaskRun> fairBuffers = new ArrayList<>()

    private int currentEmission

    private Boolean isFair0

    private TaskArrayCollector arrayCollector

    private CompilerConfiguration compilerConfig() {
        final config = new CompilerConfiguration()
        config.addCompilationCustomizers( new ASTTransformationCustomizer(TaskTemplateVarsXform) )
        config.addCompilationCustomizers( new ASTTransformationCustomizer(TaskCmdXform) )
        return config
    }

    @TestOnly
    static void reset() {
        processCount=0
        TaskRunner.errorShown.set(false)
        currentProcessor0 = null
    }

    /*
     * Initialise the process ID
     *
     * Note: processes are create in a sequential manner (by the main thread that parse the script)
     * so it does not require a synchronized block
     */
    {
        id = ++processCount
        grengine = session && session.classLoader ? new Grengine(session.classLoader, compilerConfig()) : new Grengine(compilerConfig())
        currentProcessor0 = this
    }

    @TestOnly
    protected TaskProcessor() {}

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
    TaskProcessor(String name, Executor executor, Session session, BaseScript script, ProcessConfig config, BodyDef taskBody ) {
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
        this.maxForks = config.maxForks && config.maxForks>0 ? config.maxForks as int : 0
        this.forksCount = maxForks ? new LongAdder() : null
        this.isFair0 = config.getFair()
        final arraySize = config.getArray()
        this.arrayCollector = createArrayCollector(arraySize)
        log.debug "Creating process '$name': maxForks=${maxForks}; fair=${isFair0}; array=${arraySize}"
    }

    @PackageScope TaskArrayCollector getArrayCollector() { arrayCollector }

    /**
     * @return The {@link TaskRunner} which executes the tasks of this process
     */
    TaskRunner getRunner() { runner }

    private TaskArrayCollector createArrayCollector(int arraySize) {
        if( arraySize > 0 && executor instanceof TaskArrayExecutor )
            return new TaskArrayCollector(this, executor, arraySize)
        if( arraySize > 0 )
            log.warn "Executor '${executor.name}' does not support job arrays -- the array directive will be ignored for process '$name'"
        return null
    }

    /**
     * @return The processor unique id
     */
    int getId() { id }

    /**
     * @return The {@code TaskConfig} object holding the task configuration properties
     */
    ProcessConfig getConfig() { config }

    private ProcessConfigV1 configV1() { config as ProcessConfigV1 }

    private ProcessConfigV2 configV2() { config as ProcessConfigV2 }

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

    boolean isFusionEnabled() { executor?.isFusionEnabled() ?: false }

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
    BodyDef getTaskBody() { taskBody }

    Set<String> getDeclaredNames() {
        final result = new HashSet<String>(20)
        if( config instanceof ProcessConfigV1 ) {
            result.addAll(config.getInputs().getNames())
            result.addAll(config.getOutputs().getNames())
        }
        else if( config instanceof ProcessConfigV2 ) {
            result.addAll(config.getInputs().getParams()*.getName())
            result.addAll(config.getOutputs().getParams()*.getName())
        }
        return result
    }

    LongAdder getForksCount() { forksCount }

    int getMaxForks() { maxForks }

    boolean hasErrors() { runner.getErrorCount()>0 }

    boolean isSingleton() { singleton }

    /**
     * Create a "preview" for a task run. This method is only meant for the creation of "mock" task run
     * to allow the access for the associated {@link TaskConfig} during a pipeline "preview" execution.
     *
     * Note this returns an "eventually" task configuration object. Also Inputs and output parameters are NOT
     * resolved by this method.
     *
     * @return A {@link TaskRun} object holding a reference to the associated {@link TaskConfig}
     */
    TaskRun createTaskPreview() {
        final task = new TaskRun(
                processor: this,
                type: scriptType,
                config: config.createTaskConfig(),
                context: new TaskContext(this)
        )
        task.config.context = task.context
        task.config.process = task.processor.name
        task.config.executor = task.processor.executor.name

        return task
    }

    protected void checkWarn(String msg, Map opts=null) {
        if( NF.isSyntaxParserV2() || NF.isStrictMode() )
            throw new ProcessUnrecoverableException(msg)
        if( opts )
            log.warn1(opts, msg)
        else
            log.warn(msg)
    }

    protected boolean allScalarValues

    protected boolean hasEachParams

    def run() {
        if( config instanceof ProcessConfigV1 )
            runV1()
        else if( config instanceof ProcessConfigV2 )
            runV2()
    }

    def runV1() {
        final config = configV1()

        // -- check that the task has a body
        if ( !taskBody )
            throw new IllegalStateException("Missing task body for process `$name`")

        // -- check that input tuple defines at least two elements
        def invalidInputTuple = config.getInputs().find { it instanceof TupleInParam && it.inner.size()<2 }
        if( invalidInputTuple )
            checkWarn "Input `tuple` must define at least two elements -- Check process `$name`"

        // -- check that output tuple defines at least two elements
        def invalidOutputTuple = config.getOutputs().find { it instanceof TupleOutParam && it.inner.size()<2 }
        if( invalidOutputTuple )
            checkWarn "Output `tuple` must define at least two elements -- Check process `$name`"

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
        createStateObj()

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

    protected void createStateObj() {
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
    }

    protected void createOperator() {
        def config = configV1()
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

        /**
         * The thread pool used by GPars. The thread pool to be used is set in the static
         * initializer of {@link nextflow.cli.CmdRun} class. See also {@link nextflow.util.CustomPoolFactory}
         */
        final PGroup group = Dataflow.retrieveCurrentDFPGroup()

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
            def DataflowOperator op1
            def stopAfterFirstRun = allScalarValues
            def interceptor = new BaseProcessInterceptor(opInputs, stopAfterFirstRun)
            def params = [inputs: opInputs, outputs: linkingChannels, maxForks: 1, listeners: [interceptor]]
            session.allOperators << (op1 = new DataflowOperator(group, params, forwarder))
            // fix issue #41
            start(op1)

            // set as next inputs the result channels of the iteration process
            // adding the 'control' channel removed previously
            opInputs = new ArrayList(size+1)
            opInputs.addAll( linkingChannels )
            opInputs.add( config.getInputs().getChannels().last() )
        }

        /*
         * finally create the operator
         */
        // note: do not specify the output channels in the operator declaration
        // this allows us to manage them independently from the operator life-cycle
        this.singleton = allScalarValues && !hasEachParams
        this.openPorts = createPortsArray(opInputs.size())
        config.getOutputs().setSingleton(singleton)
        def interceptor = new TaskProcessorInterceptor(opInputs, singleton)
        def params = [inputs: opInputs, maxForks: session.poolSize, listeners: [interceptor] ]
        def invoke = new InvokeTaskAdapter(this, opInputs.size())
        session.allOperators << (operator = new DataflowOperator(group, params, invoke))

        // notify the creation of a new vertex the execution DAG
        NodeMarker.addProcessNode(this, config.getInputs(), config.getOutputs())

        // fix issue #41
        start(operator)
    }

    private start(DataflowProcessor op) {
        session.addIgniter {
            log.debug "Starting process > $name"
            op.start()
        }
    }

    private AtomicIntegerArray createPortsArray(int size) {
        def result = new AtomicIntegerArray(size)
        for( int i=0; i<size; i++ )
            result.set(i, 1)
        return result
    }

    void runV2() {
        // -- check that the task has a body
        if ( !taskBody )
            throw new IllegalStateException("Missing task body for process `$name`")

        // create the state agent
        createStateObj()

        // register the processor
        session.processRegister(this)

        // determine whether the process is executed only once
        final inputs = config.getInputs().getChannels()
        this.singleton = config.getInputs().isSingleton()

        // create inputs with control channel
        final control = CH.queue()
        control.bind(Boolean.TRUE)

        final opInputs = inputs + [control]
        this.openPorts = createPortsArray(opInputs.size())

        // The thread pool used by GPars. The thread pool to be used is set in the static
        // initializer of {@link nextflow.cli.CmdRun} class. See also {@link nextflow.util.CustomPoolFactory}
        final group = Dataflow.retrieveCurrentDFPGroup()

        // note: do not specify the output channels in the operator declaration
        // this allows us to manage them independently from the operator life-cycle
        final interceptor = new TaskProcessorInterceptor(opInputs, singleton)
        final params = [inputs: opInputs, maxForks: session.poolSize, listeners: [interceptor] ]
        final invoke = new InvokeTaskAdapter(this, opInputs.size())
        this.operator = new DataflowOperator(group, params, invoke)
        session.allOperators << operator

        // notify the creation of a new vertex the execution DAG
        NodeMarker.addProcessNode(this, configV2().getInputs(), configV2().getOutputs())

        // start the operator
        start(operator)

        session.notifyProcessCreate(this)
    }

    /**
     * The processor execution body
     *
     * @param args
     *      The args array is expected to be composed by two elements:
     *      the first must be an object object of type {@link TaskStartParams},
     *      the second is the list of task input messages as received by the process
     */
    final protected void invokeTask( Object[] args ) {
        assert args.size()==2
        runner.submit((TaskStartParams) args[0], (List) args[1])
    }

    /**
     * @return A string 'she-bang' formatted to the added on top script to be executed.
     * The interpreter to be used define by the *taskConfig* property {@code shell}
     */
    static String shebangLine(shell) {
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
        result << script.stripIndent(true).trim()
        result << '\n'

        if( result[0] != '#' || result[1] != '!') {
            result.insert(0, shebangLine(shell) + '\n')
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

        if( script && script[0] == '#' && script[1] == '!') {
            return script.readLines()[0].substring(2)
        }

        return null
    }

    protected String safeTaskName(TaskRun task)  {
        return task!=null
                ? task.lazyName()
                : name
    }


    /**
     * Send a poison pill over all the outputs channel
     */
    final protected synchronized void sendPoisonPill() {
        log.trace "<$name> Sending a poison pill(s)"

        final channels = config instanceof ProcessConfigV2
            ? configV2().getOutputs().getChannels()
            : configV1().getOutputs().getChannels()

        for( DataflowWriteChannel channel : channels ){

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

    /**
     * Publish output files to a specified target folder
     *
     * @param task The task whose outputs need to be published
     * @param overwrite When {@code true} any existing file will be overwritten, otherwise the publishing is ignored
     */
    @CompileStatic
    protected void publishOutputs( TaskRun task ) {
        final publishList = task.config.getPublishDir()
        if( !publishList ) {
            return
        }

        for( PublishDir pub : publishList ) {
            publishOutputs0(task, pub)
        }
    }

    @CompileStatic
    private void publishOutputs0( TaskRun task, PublishDir publish ) {

        if( publish.overwrite == null ) {
            publish.overwrite = !task.cached
        }

        final files = getPublishFiles(task)

        publish.apply(files, task)
    }

    @CompileStatic
    private Set<Path> getPublishFiles(TaskRun task) {
        if( config instanceof ProcessConfigV1 )
            return getPublishFilesV1(task)
        if( config instanceof ProcessConfigV2 )
            return task.outputFiles
        return null
    }

    @CompileStatic
    private Set<Path> getPublishFilesV1(TaskRun task) {
        HashSet<Path> files = []
        def outputs = task.getOutputsByType(FileOutParam)
        for( Map.Entry entry : outputs ) {
            final value = entry.value
            if( value instanceof Path ) {
                files.add((Path)value)
            }
            else if( value instanceof Collection<Path> ) {
                files.addAll(value)
            }
            else if( value != null ) {
                throw new IllegalArgumentException("Unknown output file object [${value.class.name}]: ${value}")
            }
        }
        return files
    }

    /**
     * Bind the expected output files to the corresponding output channels
     * @param processor
     */
    synchronized protected void bindOutputs( TaskRun task ) {

        // bind the output
        if( isFair0 ) {
            fairBindOutputs0(task)
        }
        else {
            bindOutputs0(task)
        }

        // -- finally prints out the task output when 'debug' is true
        if( task.config.debug ) {
            task.echoStdout(session)
        }
    }

    protected void fairBindOutputs0(TaskRun task) {
        synchronized (isFair0) {
            // decrement -1 because tasks are 1-based
            final index = task.index-1
            // store the task emission values in a buffer
            fairBuffers[index-currentEmission] = task
            // check if the current task index matches the expected next emission index
            if( currentEmission == index ) {
                while( task!=null ) {
                    // bind the emission values
                    bindOutputs0(task)
                    // remove the head and try with the following
                    fairBuffers.remove(0)
                    // increase the index of the next emission
                    currentEmission++
                    // take the next task
                    task = fairBuffers[0]
                }
            }
        }
    }

    protected void bindOutputs0(TaskRun task) {
        if( config instanceof ProcessConfigV2 )
            bindOutputsV2(task)
        else if( config instanceof ProcessConfigV1 )
            bindOutputsV1(task)
    }

    @CompileStatic
    protected void bindOutputsV2(TaskRun task) {
        for( final param : configV2().getOutputs().getParams() ) {
            final value = task.outputs[param]

            log.trace "Process $name > Emitting output: ${param.name} = ${value}"
            param.getChannel().bind(value)
        }
    }

    protected void bindOutputsV1(TaskRun task) {

        // -- creates the map of all tuple values to bind
        Map<Short,List> tuples = [:]
        for( OutParam param : configV1().getOutputs() ) {
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

            case EnvOutParam:
            case ValueOutParam:
            case DefaultOutParam:
                log.trace "Process $name > collecting out param: ${param} = $value"
                tuples[param.index].add(value)
                break

            default:
                throw new IllegalArgumentException("Illegal output parameter type: $param")
            }
        }

        // -- bind out the collected values
        for( OutParam param : configV1().getOutputs() ) {
            final outValue = tuples[param.index]
            if( outValue == null )
                throw new IllegalStateException()

            if( outValue instanceof MissingParam ) {
                log.debug "Process $name > Skipping output binding because one or more optional files are missing: $outValue.missing"
                continue
            }

            log.trace "Process $name > Binding out param: ${param} = ${outValue}"
            bindOutParam(param, outValue)
        }
    }

    protected void bindOutParam( OutParam param, List values ) {
        log.trace "<$name> Binding param $param with $values"
        final x = values.size() == 1 ? values[0] : values
        final ch = param.getOutChannel()
        if( ch != null ) {
            // create a copy of the output list of operation made by a downstream task
            // can modify the list which is used internally by the task processor
            // and result in a potential error. See https://github.com/nextflow-io/nextflow/issues/3768
            final copy = x instanceof List && x instanceof Cloneable ? x.clone() : x
            // emit the final value
            ch.bind(copy)
        }
    }


    @Memoized
    ResourcesBundle getModuleBundle() {
        final script = this.getOwnerScript()
        final meta = ScriptMeta.get(script)
        // No script meta registered (e.g. processors not tied to a loaded script): nothing to resolve.
        if( meta == null )
            return null
        // Resolve the bundle when the owner script is either an included module,
        // or the entry script of a `nextflow module run` invocation (see #7087).
        return (meta.isModule() || session.isModuleRun()) ? meta.getModuleBundle() : null
    }

    @Memoized
    protected List<Path> getBinDirs() {
        final result = new ArrayList(10)
        // module bundle bin dir have priority, add before
        final bundle = session.enableModuleBinaries() ? getModuleBundle() : null
        if( bundle!=null )
            result.addAll(bundle.getBinDirs())
        // then add project bin dir
        if( executor.binDir )
            result.add(executor.binDir)
        return result
    }

    @Memoized
    boolean isLocalWorkDir() {
        return executor.workDir.fileSystem == FileSystems.default
    }

    /**
     * @return The map holding the shell environment variables for the task to be executed
     */
    @Memoized
    Map<String,String> getProcessEnvironment() {

        def result = new LinkedHashMap<String,String>(20)

        // add the taskConfig environment entries
        if( session.config.env instanceof Map ) {
            session.config.env.each { name, value ->
                // skip entries with a null value (e.g. `env(HOST_VAR)` referencing a
                // host environment variable that is not set) instead of exporting
                // them as an empty string with a spurious warning -- see #5722
                if( value != null )
                    result.put( name, value.toString() )
            }
        }
        else {
            log.debug "Invalid 'session.config.env' object: ${session.config.env?.class?.name}"
        }

        // append the 'bin' folder to the task environment
        List<Path> paths
        if( isLocalWorkDir() && (paths=getBinDirs()) ) {
            for( Path it : paths ) {
                if( result.containsKey('PATH') ) {
                    // note: do not escape potential blanks in the bin path because the PATH
                    // variable is enclosed in `"` when in rendered in the launcher script -- see #630
                    result['PATH'] =  "${result['PATH']}:${it}".toString()
                }
                else {
                    // note: append custom bin path *after* the system PATH
                    // to prevent unnecessary network round-trip for each command
                    // when the added path is a shared file system directory
                    result['PATH'] = "\$PATH:${it}".toString()
                }
            }
        }

        return Collections.unmodifiableMap(result)
    }

    /**
     * Given a map holding variables key-value pairs, create a script fragment
     * exporting the required environment variables
     */
    static String bashEnvironmentScript( Map<String,String> environment, boolean escape=false ) {
        if( !environment )
            return null

        final List script = []
        for( String name : environment.keySet() ) {
            String value = environment.get(name)
            if( !ENV_VAR_NAME.matcher(name).matches() )
                log.trace "Illegal environment variable name: '${name}' -- This variable definition is ignored"
            else if( !value ) {
                log.warn "Environment variable `$name` evaluates to an empty value"
                script << "export $name=''"
            }
            else if( !escape ) {
                script << /export $name="$value"/
            }
            else {
                // escape both wrapping double quotes and the dollar var placeholder
                script << /export $name="${Escape.variable(value)}"/
            }
        }
        script << ''

        return script.join('\n')
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
     * Invoked when a task execution has completed, that is when the task outputs
     * have been resolved, or the task has failed, or its execution was skipped
     * because of a `when` guard, a store directory or a cached result.
     *
     * This is the boundary between the execution of a task and the dataflow
     * interface of the process: it is the only point at which task results are
     * emitted to the process output channels.
     *
     * @param task The {@code TaskRun} instance that has completed
     */
    protected void taskCompleted( TaskRun task ) {
        log.trace "Finalize process > ${safeTaskName(task)}"

        // -- bind output (files)
        if( task.canBind ) {
            bindOutputs(task)
            publishOutputs(task)
        }

        // increment the number of processes executed
        state.update { StateObj it -> it.incCompleted() }
    }

    protected void closeProcess() {
        arrayCollector?.close()
    }

    protected void terminateProcess() {
        log.trace "<${name}> Sending poison pills and terminating process"
        sendPoisonPill()
        session.notifyProcessTerminate(this)
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
        if( config instanceof ProcessConfigV1 )
            dumpOpenPorts(result)

        return result.toString()
    }

    void dumpOpenPorts(StringBuilder result) {
        for( int i=0; i<openPorts.length(); i++ ) {
            def last = i == openPorts.length()-1
            def param = configV1().getInputs()[i]
            def chnnl = param?.inChannel
            def isValue = chnnl instanceof DataflowExpression
            def type = last ? '(cntrl)' : (isValue ? '(value)' : '(queue)')
            def channel = param && !(param instanceof TupleInParam) ? param.getName() : '-'
            def status; if( isValue ) { status = !chnnl.isBound() ? 'OPEN  ' : 'bound ' }
            else status = type == '(queue)' ? (openPorts.get(i) ? 'OPEN  ' : 'closed') : '-     '
            result << "  port $i: $type ${status}; channel: $channel\n"
        }
    }

    /*
     * logger class for the *iterator* processor
     */
    class BaseProcessInterceptor extends DataflowEventAdapter {

        final List<DataflowReadChannel> inputs

        final boolean stopAfterFirstRun

        final int len

        final DataflowQueue control

        final int first

        BaseProcessInterceptor( List<DataflowReadChannel> inputs, boolean stop ) {
            this.inputs = new ArrayList<>(inputs)
            this.stopAfterFirstRun = stop
            this.len = inputs.size()
            this.control = (DataflowQueue)inputs.get(len-1)
            this.first = inputs.findIndexOf { CH.isChannelQueue(it) }
        }

        @Override
        Object messageArrived(final DataflowProcessor processor, final DataflowReadChannel<Object> channel, final int index, final Object message) {
            if( len == 1 || stopAfterFirstRun ) {
                // -- kill itself
                control.bind(PoisonPill.instance)
            }
            else if( index == first ) {
                // the `if` condition guarantees only and only one signal message (the true value)
                // is bound to the control message for a complete set of input values delivered
                // to the process -- the control message is need to keep the process running
                control.bind(Boolean.TRUE)
            }

            return message
        }
    }

    /**
     *  Intercept dataflow process events
     */
    class TaskProcessorInterceptor extends BaseProcessInterceptor {

        TaskProcessorInterceptor(List<DataflowReadChannel> inputs, boolean stop) {
            super(inputs, stop)
        }

        @Override
        List<Object> beforeRun(final DataflowProcessor processor, final List<Object> messages) {
            // apparently auto if-guard instrumented by @Slf4j is not honoured in inner classes - add it explicitly
            if( log.isTraceEnabled() )
                log.trace "<${name}> Before run -- messages: ${messages}"
            // the counter must be incremented here, otherwise it won't be consistent
            state.update { StateObj it -> it.incSubmitted() }
            // task index must be created here to guarantee consistent ordering
            // with the sequence of messages arrival since this method is executed in a thread safe manner
            final params = new TaskStartParams(TaskId.next(), indexCount.incrementAndGet())
            final result = new ArrayList(2)
            result[0] = params
            result[1] = messages
            return result
        }


        @Override
        void afterRun(DataflowProcessor processor, List<Object> messages) {
            // apparently auto if-guard instrumented by @Slf4j is not honoured in inner classes - add it explicitly
            if( log.isTraceEnabled() )
                log.trace "<${name}> After run"
            currentTask.remove()
        }

        @Override
        Object messageArrived(final DataflowProcessor processor, final DataflowReadChannel<Object> channel, final int index, final Object message) {
            // apparently auto if-guard instrumented by @Slf4j is not honoured in inner classes - add it explicitly
            if( log.isTraceEnabled() && config instanceof ProcessConfigV1 ) {
                def channelName = configV1().getInputs()?.names?.get(index)
                def taskName = currentTask.get()?.name ?: name
                log.trace "<${taskName}> Message arrived -- ${channelName} => ${message}"
            }

            super.messageArrived(processor, channel, index, message)
        }

        @Override
        Object controlMessageArrived(final DataflowProcessor processor, final DataflowReadChannel<Object> channel, final int index, final Object message) {
            // apparently auto if-guard instrumented by @Slf4j is not honoured in inner classes - add it explicitly
            if( log.isTraceEnabled() && config instanceof ProcessConfigV1 ) {
                def channelName = configV1().getInputs()?.names?.get(index)
                def taskName = currentTask.get()?.name ?: name
                log.trace "<${taskName}> Control message arrived ${channelName} => ${message}"
            }

            super.controlMessageArrived(processor, channel, index, message)

            if( message == PoisonPill.instance ) {
                // apparently auto if-guard instrumented by @Slf4j is not honoured in inner classes - add it explicitly
                if( log.isTraceEnabled() )
                    log.trace "<${name}> Poison pill arrived; port: $index"
                openPorts.set(index, 0) // mark the port as closed
                state.update { StateObj it -> it.poison() }
            }

            return message
        }

        @Override
        void afterStop(final DataflowProcessor processor) {
            // apparently auto if-guard instrumented by @Slf4j is not honoured in inner classes - add it explicitly
            if( log.isTraceEnabled() )
                log.trace "<${name}> After stop"
            closeProcess()
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
        boolean onException(final DataflowProcessor processor, final Throwable error) {
            // return `true` to terminate the dataflow processor
            runner.handleException( error, currentTask.get() )
        }
    }
}
