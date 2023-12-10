/*
 * Copyright 2013-2023, Seqera Labs
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

import static nextflow.processor.ErrorStrategy.*

import java.lang.reflect.InvocationTargetException
import java.nio.file.FileSystems
import java.nio.file.NoSuchFileException
import java.nio.file.Path
import java.nio.file.Paths
import java.util.concurrent.atomic.AtomicBoolean
import java.util.concurrent.atomic.AtomicInteger
import java.util.concurrent.atomic.LongAdder
import java.util.regex.Matcher
import java.util.regex.Pattern

import ch.artecat.grengine.Grengine
import com.google.common.hash.HashCode
import groovy.json.JsonOutput
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
import nextflow.Nextflow
import nextflow.Session
import nextflow.ast.DslCodeVisitor
import nextflow.ast.TaskCmdXform
import nextflow.ast.TaskTemplateVarsXform
import nextflow.cloud.CloudSpotTerminationException
import nextflow.dag.NodeMarker
import nextflow.exception.FailedGuardException
import nextflow.exception.IllegalArityException
import nextflow.exception.MissingFileException
import nextflow.exception.MissingValueException
import nextflow.exception.ProcessException
import nextflow.exception.ProcessFailedException
import nextflow.exception.ProcessRetryableException
import nextflow.exception.ProcessSubmitTimeoutException
import nextflow.exception.ProcessUnrecoverableException
import nextflow.exception.ShowOnlyExceptionMessage
import nextflow.exception.UnexpectedException
import nextflow.executor.CachedTaskHandler
import nextflow.executor.Executor
import nextflow.executor.StoredTaskHandler
import nextflow.extension.CH
import nextflow.extension.DataflowHelper
import nextflow.file.FileHelper
import nextflow.file.FileHolder
import nextflow.file.FilePorter
import nextflow.script.BaseScript
import nextflow.script.BodyDef
import nextflow.script.ProcessConfig
import nextflow.script.ScriptMeta
import nextflow.script.ScriptType
import nextflow.script.TaskClosure
import nextflow.script.bundle.ResourcesBundle
import nextflow.util.ArrayBag
import nextflow.util.BlankSeparatedList
import nextflow.util.CacheHelper
import nextflow.util.Escape
import nextflow.util.LazyHelper
import nextflow.util.LockManager
import nextflow.util.LoggerHelper
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

    final private static Pattern QUESTION_MARK = ~/(\?+)/

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
     * Count the number of time an error occurred
     */
    private volatile int errorCount


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
     * Whenever the process is closed (ie. received the STOP signal)
     */
    protected AtomicBoolean closed = new AtomicBoolean(false)

    /**
     * Process ID number. The first is 1, the second 2 and so on ..
     */
    private final int id

    private LongAdder forksCount

    private int maxForks

    private static int processCount

    private static LockManager lockManager = new LockManager()

    private List<Map<Short,List>> fairBuffers = new ArrayList<>()

    private int currentEmission

    private Boolean isFair0

    private CompilerConfiguration compilerConfig() {
        final config = new CompilerConfiguration()
        config.addCompilationCustomizers( new ASTTransformationCustomizer(TaskTemplateVarsXform) )
        config.addCompilationCustomizers( new ASTTransformationCustomizer(TaskCmdXform) )
        return config
    }

    @TestOnly
    static void reset() {
        processCount=0
        errorShown.set(false)
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
        this.maxForks = config.maxForks ? config.maxForks as int : 0
        this.forksCount = maxForks ? new LongAdder() : null
        this.isFair0 = config.getFair()
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
    BodyDef getTaskBody() { taskBody }

    Set<String> getDeclaredNames() {
        Set<String> result = new HashSet<>(20)
        result.addAll(config.getInputs().getNames())
        result.addAll(config.getOutputs().getNames())
        return result
    }

    LongAdder getForksCount() { forksCount }

    int getMaxForks() { maxForks }

    boolean hasErrors() { errorCount>0 }

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
        if( NF.isStrictMode() )
            throw new ProcessUnrecoverableException(msg)
        if( opts )
            log.warn1(opts, msg)
        else
            log.warn(msg)
    }

    def run(DataflowReadChannel source) {

        // -- check that the task has a body
        if ( !taskBody )
            throw new IllegalStateException("Missing task body for process `$name`")

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
        createOperator(source)

        session.notifyProcessCreate(this)

        /*
         * When there is a single output channel, return let returns that item
         * otherwise return the list
         */
        final result = config.getOutputs().getChannels()
        return result.size() == 1 ? result[0] : result
    }

    protected void createOperator(DataflowReadChannel source) {
        // determine whether the process is executed only once
        this.singleton = !CH.isChannelQueue(source)

        // create inputs with control channel
        final control = CH.queue()
        control.bind(Boolean.TRUE)

        final opInputs = [source, control]

        /**
         * The thread pool used by GPars. The thread pool to be used is set in the static
         * initializer of {@link nextflow.cli.CmdRun} class. See also {@link nextflow.util.CustomPoolFactory}
         */
        final PGroup group = Dataflow.retrieveCurrentDFPGroup()

        // note: do not specify the output channels in the operator declaration
        // this allows us to manage them independently from the operator life-cycle
        def interceptor = new TaskProcessorInterceptor(source, control, singleton)
        def params = [inputs: opInputs, maxForks: session.poolSize, listeners: [interceptor] ]
        this.operator = new DataflowOperator(group, params, this.&invokeTask)

        // notify the creation of a new vertex the execution DAG
        NodeMarker.addProcessNode(this, config.getInputs(), config.getOutputs())

        // start the operator
        session.allOperators << operator
        session.addIgniter {
            log.debug "Starting process > $name"
            operator.start()
        }
    }

    final protected void invokeTask( TaskStartParams params, List values ) {
        // create and initialize the task instance to be executed
        log.trace "Invoking task > $name with params=$params; values=$values"

        // -- create the task run instance
        final task = createTaskRun(params)
        // -- set the task instance as the current in this thread
        currentTask.set(task)

        // -- map the inputs to a map and use to delegate closure values interpolation
        resolveTaskInputs(task, values)

        // verify that `when` guard, when specified, is satisfied
        if( !checkWhenGuard(task) )
            return

        TaskClosure block
        if( session.stubRun && (block=task.config.getStubBlock()) ) {
            task.resolve(block)
        }
        else {
            // -- prepend task config to arguments (for process function)
            values.push(task.config)

            // -- resolve the task command script
            task.resolve(taskBody, values.toArray())
        }

        // -- verify if exists a stored result for this case,
        //    if true skip the execution and return the stored data
        if( checkStoredOutput(task) )
            return

        def hash = createTaskHashKey(task)
        checkCachedOrLaunchTask(task, hash, resumable)
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

    final protected TaskRun createTaskRun(TaskStartParams params) {
        final task = new TaskRun(
                id: params.id,
                index: params.index,
                processor: this,
                type: scriptType,
                config: config.createTaskConfig(),
                context: new TaskContext(this)
        )

        // setup config
        task.config.index = task.index
        task.config.process = task.processor.name
        task.config.executor = task.processor.executor.name

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
    @CompileStatic
    final protected void checkCachedOrLaunchTask( TaskRun task, HashCode hash, boolean shouldTryCache ) {

        int tries = task.failCount +1
        while( true ) {
            hash = CacheHelper.defaultHasher().newHasher().putBytes(hash.asBytes()).putInt(tries).hash()

            Path resumeDir = null
            boolean exists = false
            try {
                final entry = session.cache.getTaskEntry(hash, this)
                resumeDir = entry ? FileHelper.asPath(entry.trace.getWorkDir()) : null
                if( resumeDir )
                    exists = resumeDir.exists()

                log.trace "[${safeTaskName(task)}] Cacheable folder=${resumeDir?.toUriString()} -- exists=$exists; try=$tries; shouldTryCache=$shouldTryCache; entry=$entry"
                final cached = shouldTryCache && exists && entry.trace.isCompleted() && checkCachedOutput(task.clone(), resumeDir, hash, entry)
                if( cached )
                    break
            }
            catch (Throwable t) {
                log.warn1("[${safeTaskName(task)}] Unable to resume cached task -- See log file for details", causedBy: t)
            }

            if( exists ) {
                tries++
                continue
            }

            final lock = lockManager.acquire(hash)
            final workDir = task.getWorkDirFor(hash)
            try {
                if( resumeDir != workDir )
                    exists = workDir.exists()
                if( !exists && !workDir.mkdirs() )
                    throw new IOException("Unable to create directory=$workDir -- check file system permissions")
            }
            finally {
                lock.release()
            }

            // submit task for execution
            submitTask( task, hash, workDir )
            break
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
            log.trace "[${safeTaskName(task)}] Store dir not set -- return false"
            return false
        }

        if( !task.config.getStoreDir().exists() ) {
            log.trace "[${safeTaskName(task)}] Store dir does not exists > ${task.config.storeDir} -- return false"
            // no folder -> no cached result
            return false
        }


        try {
            // -- expose task exit status to make accessible as output value
            task.config.exitStatus = TaskConfig.EXIT_ZERO
            // -- check if all output resources are available
            collectOutputs(task)
            log.info "[skipping] Stored process > ${safeTaskName(task)}"
            // set the exit code in to the task object
            task.exitStatus = TaskConfig.EXIT_ZERO
            task.cached = true
            session.notifyTaskCached(new StoredTaskHandler(task))

            // -- now bind the results
            finalizeTask0(task)
            return true
        }
        catch( MissingFileException | MissingValueException e ) {
            log.trace "[${safeTaskName(task)}] Missed store > ${e.getMessage()} -- folder: ${task.config.storeDir}"
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
    final boolean checkCachedOutput(TaskRun task, Path folder, HashCode hash, TaskEntry entry) {

        // check if exists the task exit code file
        def exitCode = null
        def exitFile = folder.resolve(TaskRun.CMD_EXIT)
        if( task.type == ScriptType.SCRIPTLET ) {
            def str
            try {
                str = exitFile.text?.trim()
            }
            catch( IOException e ) {
                log.trace "[${safeTaskName(task)}] Exit file can't be read > $exitFile -- return false -- Cause: ${e.message}"
                return false
            }

            exitCode = str.isInteger() ? str.toInteger() : null
            if( !task.isSuccess(exitCode) ) {
                log.trace "[${safeTaskName(task)}] Exit code is not valid > $str -- return false"
                return false
            }
        }

        /*
         * verify cached context map
         */
        if( !entry ) {
            log.trace "[${safeTaskName(task)}] Missing cache entry -- return false"
            return false
        }

        if( !entry.context ) {
            log.trace "[${safeTaskName(task)}] Missing cache context -- return false"
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
            // -- set task properties in order to resolve outputs
            task.workDir = folder
            task.stdout = stdoutFile
            task.config.exitStatus = exitCode
            // -- check if all output resources are available
            collectOutputs(task)

            // set the exit code in to the task object
            task.cached = true
            task.hash = hash
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
            log.trace "[${safeTaskName(task)}] Missed cache > ${e.getMessage()} -- folder: $folder"
            task.exitStatus = Integer.MAX_VALUE
            task.workDir = null
            task.stdout = null
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
        log.trace "Task fault (2): $fault"

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
        log.debug "Handling unexpected condition for\n  task: name=${safeTaskName(task)}; work-dir=${task?.workDirStr}\n  error [${error?.class?.name}]: ${error?.getMessage()?:error}"

        ErrorStrategy errorStrategy = TERMINATE
        final List<String> message = []
        try {
            // -- do not recoverable error, just re-throw it
            if( error instanceof Error ) throw error

            // -- retry without increasing the error counts
            if( task && (error.cause instanceof ProcessRetryableException || error.cause instanceof CloudSpotTerminationException) ) {
                if( error.cause instanceof ProcessRetryableException )
                    log.info "[$task.hashLog] NOTE: ${error.message} -- Execution is retried"
                else
                    log.info "[$task.hashLog] NOTE: ${error.message} -- Cause: ${error.cause.message} -- Execution is retried"
                task.failCount+=1
                final taskCopy = task.makeCopy()
                session.getExecService().submit {
                    try {
                        taskCopy.runType = RunType.RETRY
                        checkCachedOrLaunchTask( taskCopy, taskCopy.hash, false )
                    }
                    catch( Throwable e ) {
                        log.error("Unable to re-submit task `${taskCopy.name}`", e)
                        session.abort(e)
                    }
                }
                task.failed = true
                task.errorAction = RETRY
                return RETRY
            }

            final submitTimeout = error.cause instanceof ProcessSubmitTimeoutException
            final submitErrMsg = submitTimeout ? error.cause.message : null
            final int submitRetries = submitTimeout ? ++task.submitRetries : 0
            final int taskErrCount = !submitTimeout && task ? ++task.failCount : 0
            final int procErrCount = !submitTimeout ? ++errorCount : errorCount

            // -- when is a task level error and the user has chosen to ignore error,
            //    just report and error message and DO NOT stop the execution
            if( task && error instanceof ProcessException ) {
                // expose current task exit status
                task.config.exitStatus = task.exitStatus
                task.config.errorCount = procErrCount
                task.config.retryCount = taskErrCount

                errorStrategy = checkErrorStrategy(task, error, taskErrCount, procErrCount, submitRetries)
                if( errorStrategy.soft ) {
                    def msg = "[$task.hashLog] NOTE: ${submitTimeout ? submitErrMsg : error.message}"
                    if( errorStrategy == IGNORE ) msg += " -- Error is ignored"
                    else if( errorStrategy == RETRY ) msg += " -- Execution is retried (${submitTimeout ? submitRetries : taskErrCount})"
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
                log.trace "Task errorShown=${errorShown.get()}; aborted=${session.aborted}"
                return errorStrategy
            }

            def dumpStackTrace = log.isTraceEnabled()
            message << "Error executing process > '${safeTaskName(task)}'"
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

    protected String safeTaskName(TaskRun task)  {
        return task!=null
                ? task.lazyName()
                : name
    }

    protected ErrorStrategy checkErrorStrategy( TaskRun task, ProcessException error, final int taskErrCount, final int procErrCount, final submitRetries ) {

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

            if( (procErrCount < maxErrors || maxErrors == -1) && taskErrCount <= maxRetries && submitRetries <= maxRetries ) {
                final taskCopy = task.makeCopy()
                session.getExecService().submit({
                    try {
                        taskCopy.config.attempt = taskErrCount+1
                        taskCopy.config.submitAttempt = submitRetries+1
                        taskCopy.runType = RunType.RETRY
                        taskCopy.resolve(taskBody)
                        checkCachedOrLaunchTask( taskCopy, taskCopy.hash, false )
                    }
                    catch( Throwable e ) {
                        log.error("Unable to re-submit task `${safeTaskName(taskCopy)}`", e)
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
            error.source.stripIndent(true).eachLine {
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
            task.script?.stripIndent(true)?.trim()?.eachLine {
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
            for( String it : lines ) {
                message << "  ${stripWorkDir(it, task.workDir)}"
            }

            // - the tail of the process stderr
            lines = task.dumpStderr(max)
            if( lines ) {
                message << "\nCommand error:"
                for( String it : lines ) {
                    message << "  ${stripWorkDir(it, task.workDir)}"
                }
            }
            // - this is likely a task wrapper issue
            else if( task.exitStatus != 0 ) {
                lines = task.dumpLogFile(max)
                if( lines ) {
                    message << "\nCommand wrapper:"
                    for( String it : lines ) {
                        message << "  ${stripWorkDir(it, task.workDir)}"
                    }
                }
            }

        }
        else {
            if( task?.source )  {
                message << "Source block:"
                task.source.stripIndent(true).eachLine {
                    message << "  $it"
                }
            }

        }

        if( task?.workDir )
            message << "\nWork dir:\n  ${task.workDirStr}"

        message << "\nTip: ${getRndTip()}"

        return message
    }

    private static String stripWorkDir(String line, Path workDir) {
        if( workDir==null ) return line
        if( workDir.fileSystem != FileSystems.default ) return line
        return workDir ? line.replace(workDir.toString()+'/','') : line
    }

    static List tips = [
            'when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line',
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

        for( DataflowWriteChannel channel : config.getOutputs().getChannels() ){

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
        if( error instanceof ShowOnlyExceptionMessage || !error.cause )
            message = err0(error)
        else
            message = err0(error.cause)

        result
            .append('  ')
            .append(message)
            .append('\n')
            .toString()
    }


    static String err0(Throwable e) {
        final fail = e instanceof InvocationTargetException ? e.targetException : e

        if( fail instanceof NoSuchFileException ) {
            return "No such file or directory: $fail.message"
        }
        if( fail instanceof MissingPropertyException ) {
            def name = fail.property ?: LoggerHelper.getDetailMessage(fail)
            def result = "No such variable: ${name}"
            def details = LoggerHelper.findErrorLine(fail)
            if( details )
                result += " -- Check script '${details[0]}' at line: ${details[1]}"
            return result
        }
        def result = fail.message ?: fail.toString()
        def details = LoggerHelper.findErrorLine(fail)
        if( details ){
            result += " -- Check script '${details[0]}' at line: ${details[1]}"
        }
        return result
    }

    /**
     * Publish output files to a specified target folder
     *
     * @param task The task whose outputs need to be published
     * @param overwrite When {@code true} any existing file will be overwritten, otherwise the publishing is ignored
     */
    @CompileStatic
    protected void publishOutputs( TaskRun task ) {
        final publishers = task.config.getPublishDir()

        for( PublishDir publisher : publishers ) {
            if( publisher.overwrite == null )
                publisher.overwrite = !task.cached

            publisher.apply(task)
        }
    }

    /**
     * Bind the expected output files to the corresponding output channels
     * @param processor
     */
    synchronized protected void bindOutputs( TaskRun task ) {

        // bind the output
        if( isFair0 ) {
            fairBindOutputs0(task.outputs, task)
        }
        else {
            bindOutputs0(task.outputs)
        }

        // -- finally prints out the task output when 'debug' is true
        if( task.config.debug ) {
            task.echoStdout(session)
        }
    }

    protected void fairBindOutputs0(List emissions, TaskRun task) {
        synchronized (isFair0) {
            // decrement -1 because tasks are 1-based
            final index = task.index-1
            // store the task emission values in a buffer
            fairBuffers[index-currentEmission] = emissions
            // check if the current task index matches the expected next emission index
            if( currentEmission == index ) {
                while( emissions!=null ) {
                    // bind the emission values
                    bindOutputs0(emissions)
                    // remove the head and try with the following
                    fairBuffers.remove(0)
                    // increase the index of the next emission
                    currentEmission++
                    // take the next emissions 
                    emissions = fairBuffers[0]
                }
            }
        }
    }

    protected void bindOutputs0(List outputs) {
        // -- bind out the collected values
        for( int i = 0; i < config.getOutputs().size(); i++ ) {
            final param = config.getOutputs()[i]
            final value = outputs[i]

            if( value == null ) {
                log.debug "Process $name > Skipping output binding because one or more optional files are missing: ${param.name}"
                continue
            }

            // clone collection values before emitting them so that the task processor
            // can iterate over them without causing a race condition
            // see https://github.com/nextflow-io/nextflow/issues/3768
            log.trace "Process $name > Emitting output: ${param.name} = ${value}"
            final copy = value instanceof Collection && value instanceof Cloneable ? value.clone() : value
            param.getChannel().bind(copy)
        }
    }

    /**
     * Once the task has completed this method is invoked to collected all the task results
     *
     * @param task
     */
    protected void collectOutputs( TaskRun task ) {
        task.outputs = config.getOutputs().collect( param -> param.resolve(task) )
        task.canBind = true
    }

    @Memoized
    ResourcesBundle getModuleBundle() {
        final script = this.getOwnerScript()
        final meta = ScriptMeta.get(script)
        return meta?.isModule() ? meta.getModuleBundle() : null
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
                result.put( name, value?.toString() )
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

    protected Path normalizeToPath( obj ) {
        if( obj instanceof Path )
            return obj

        if( obj == null )
            throw new ProcessUnrecoverableException("Path value cannot be null")
        
        if( !(obj instanceof CharSequence) )
            throw new ProcessUnrecoverableException("Not a valid path value type: ${obj.getClass().getName()} ($obj)")

        def str = obj.toString().trim()
        if( str.contains('\n') )
            throw new ProcessUnrecoverableException("Path value cannot contain a new-line character: $str")
        if( str.startsWith('/') )
            return FileHelper.asPath(str)
        if( FileHelper.getUrlProtocol(str) )
            return FileHelper.asPath(str)
        if( !str )
            throw new ProcessUnrecoverableException("Path value cannot be empty")
        
        throw new ProcessUnrecoverableException("Not a valid path value: '$str'")
    }

    protected List<FileHolder> normalizeInputToFiles( Object obj, int count, boolean coerceToPath, FilePorter.Batch batch ) {

        Collection allItems = obj instanceof Collection ? obj : [obj]
        def len = allItems.size()

        // use a bag so that cache hash key is not affected by file entries order
        def files = new ArrayBag<FileHolder>(len)
        for( def item : allItems ) {

            if( item instanceof Path || coerceToPath ) {
                def path = normalizeToPath(item)
                def target = executor.isForeignFile(path) ? batch.addToForeign(path) : path
                def holder = new FileHolder(target)
                files << holder
            }
            else {
                files << normalizeInputToFile(item, "input.${++count}")
            }
        }

        return files
    }

    protected singleItemOrList( List<FileHolder> items, boolean single, ScriptType type ) {
        assert items != null

        if( items.size() == 1 && single ) {
            return makePath(items[0],type)
        }

        def result = new ArrayList(items.size())
        for( int i=0; i<items.size(); i++ ) {
            result.add( makePath(items[i],type) )
        }
        return new BlankSeparatedList(result)
    }

    private Path makePath( FileHolder holder, ScriptType type ) {
        if( type == ScriptType.SCRIPTLET ) {
            return new TaskPath(holder)
        }
        if( type == ScriptType.GROOVY) {
            // the real path for the native task needs to be fixed -- see #378
            return Paths.get(holder.stageName)
        }
        throw new IllegalStateException("Unknown task type: $type")
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
             * When name do not contain any wildcards *BUT* multiple files are provide
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

    final protected void resolveTaskInputs( TaskRun task, List values ) {

        final inputs = config.getInputs()
        final ctx = task.context

        // -- add input params to task context
        for( int i = 0; i < inputs.size(); i++ )
            ctx.put(inputs[i].getName(), values[i])

        // -- resolve local variables
        for( def entry : inputs.getVariables() )
            ctx.put(entry.key, LazyHelper.resolve(ctx, entry.value))

        // -- resolve environment vars
        for( def entry : inputs.getEnv() )
            task.env.put(entry.key, LazyHelper.resolve(ctx, entry.value))

        // -- resolve stdin
        task.stdin = LazyHelper.resolve(ctx, inputs.stdin)

        // -- resolve input files
        final allNames = new HashMap<String,Integer>()
        int count = 0
        final batch = session.filePorter.newBatch(executor.getStageDir())

        for( def param : config.getInputs().getFiles() ) {
            final val = param.resolve(ctx)
            final normalized = normalizeInputToFiles(val, count, param.isPathQualifier(), batch)
            final resolved = expandWildcards( param.getFilePattern(ctx), normalized )

            if( !param.isValidArity(resolved.size()) )
                throw new IllegalArityException("Incorrect number of input files for process `${safeTaskName(task)}` -- expected ${param.arity}, found ${resolved.size()}")

            // add to context if the path was declared with a variable name
            if( param.name )
                ctx.put( param.name, singleItemOrList(resolved, param.isSingle(), task.type) )

            count += resolved.size()

            for( FileHolder item : resolved ) {
                Integer num = allNames.getOrCreate(item.stageName, 0) +1
                allNames.put(item.stageName,num)
            }

            task.inputFiles.addAll(resolved)
        }

        // -- set the delegate map as context in the task config
        //    so that lazy directives will be resolved against it
        task.config.context = ctx

        // -- check conflicting file names
        def conflicts = allNames.findAll { name, num -> num>1 }
        if( conflicts ) {
            log.debug("Process $name > collision check staging file names: $allNames")
            def message = "Process `$name` input file name collision -- There are multiple input files for each of the following file names: ${conflicts.keySet().join(', ')}"
            throw new ProcessUnrecoverableException(message)
        }

        // -- download foreign files
        session.filePorter.transfer(batch)
    }

    final protected HashCode createTaskHashKey(TaskRun task) {

        List keys = [ session.uniqueId, name, task.source ]

        if( task.isContainerEnabled() )
            keys << task.getContainerFingerprint()

        // add task inputs
        final inputVars = getTaskInputVars(task)
        if( inputVars )
            keys.add(inputVars.entrySet())
        if( task.env )
            keys.add(task.env.entrySet())
        if( task.inputFiles )
            keys.add(task.inputFiles)
        if( task.stdin )
            keys.add(task.stdin)

        // add all variable references in the task script but not declared as input/output
        def vars = getTaskGlobalVars(task)
        if( vars ) {
            log.trace "Task: $name > Adding script vars hash code: ${vars}"
            keys.add(vars.entrySet())
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
        
        final conda = task.getCondaEnv()
        if( conda ) {
            keys.add(conda)
        }

        final spack = task.getSpackEnv()
        final arch = task.getConfig().getArchitecture()

        if( spack ) {
            keys.add(spack)

            if( arch ) {
                keys.add(arch)
            }
        }

        if( session.stubRun ) {
            keys.add('stub-run')
        }

        final mode = config.getHashMode()
        final hash = computeHash(keys, mode)
        if( session.dumpHashes ) {
            session.dumpHashes=='json'
                ? traceInputsHashesJson(task, keys, mode, hash)
                : traceInputsHashes(task, keys, mode, hash)
        }
        return hash
    }

    HashCode computeHash(List keys, CacheHelper.HashMode mode) {
        try {
            return CacheHelper.hasher(keys, mode).hash()
        }
        catch (Throwable e) {
            final msg = "Something went wrong while creating task '$name' unique id -- Offending keys: ${ keys.collect {"\n - type=${it.getClass().getName()} value=$it"} }"
            throw new UnexpectedException(msg,e)
        }
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
        List<Path> result = []
        def tokenizer = new StringTokenizer(script," \t\n\r\f()[]{};&|<>`")
        while( tokenizer.hasMoreTokens() ) {
            def token = tokenizer.nextToken()
            def path = session.binEntries.get(token)
            if( path )
                result.add(path)
        }
        return result
    }

    private void traceInputsHashesJson( TaskRun task, List entries, CacheHelper.HashMode mode, hash ) {
        final collector = (item) -> [
            hash: CacheHelper.hasher(item, mode).hash().toString(),
            type: item?.getClass()?.getName(),
            value: item?.toString()
        ]
        final json = JsonOutput.toJson(entries.collect(collector))
        log.info "[${safeTaskName(task)}] cache hash: ${hash}; mode: ${mode}; entries: ${JsonOutput.prettyPrint(json)}"
    }

    private void traceInputsHashes( TaskRun task, List entries, CacheHelper.HashMode mode, hash ) {

        def buffer = new StringBuilder()
        buffer.append("[${safeTaskName(task)}] cache hash: ${hash}; mode: $mode; entries: \n")
        for( Object item : entries ) {
            buffer.append( "  ${CacheHelper.hasher(item, mode).hash()} [${item?.getClass()?.getName()}] $item \n")
        }

        log.info(buffer.toString())
    }

    protected Map<String,Object> getTaskInputVars(TaskRun task) {
        final result = [:]
        final inputs = config.getInputs()
        final inputVars = inputs.getNames() - inputs.getFiles()*.getName()
        for( String var : inputVars )
            result.put(var, task.context.get(var))
        return result
    }

    protected Map<String,Object> getTaskGlobalVars(TaskRun task) {
        final result = task.getGlobalVars(ownerScript.getBinding())
        final directives = getTaskExtensionDirectiveVars(task)
        result.putAll(directives)
        return result
    }

    protected Map<String,Object> getTaskExtensionDirectiveVars(TaskRun task) {
        final variableNames = task.getVariableNames()
        final result = new HashMap(variableNames.size())
        final taskConfig = task.config
        for( String key : variableNames ) {
            if( !key.startsWith('task.ext.') ) continue
            final value = taskConfig.eval(key.substring(5))
            result.put(key, value)
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
        log.trace "[${safeTaskName(task)}] actual run folder: ${folder}"

        // set hash-code & working directory
        task.hash = hash
        task.workDir = folder
        task.config.workDir = folder
        task.config.hash = hash.toString()
        task.config.name = task.getName()

        // add the task to the collection of running tasks
        executor.submit(task)

    }

    protected boolean checkWhenGuard(TaskRun task) {

        try {
            def pass = task.config.getGuard(DslCodeVisitor.PROCESS_WHEN)
            if( pass ) {
                return true
            }

            log.trace "Task ${safeTaskName(task)} is not executed because `when` condition is not verified"
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
        log.trace "finalizing process > ${safeTaskName(task)} -- $task"

        def fault = null
        try {
            // -- verify task exit status
            if( task.error )
                throw new ProcessFailedException("Process `${safeTaskName(task)}` failed", task.error)

            if( task.type == ScriptType.SCRIPTLET ) {
                if( task.exitStatus == Integer.MAX_VALUE )
                    throw new ProcessFailedException("Process `${safeTaskName(task)}` terminated for an unknown reason -- Likely it has been terminated by the external system")

                if ( !task.isSuccess() )
                    throw new ProcessFailedException("Process `${safeTaskName(task)}` terminated with an error exit status (${task.exitStatus})")
            }

            // -- expose task exit status to make accessible as output value
            task.config.exitStatus = task.exitStatus
            // -- if it's OK collect results and finalize
            collectOutputs(task)
        }
        catch ( Throwable error ) {
            fault = resumeOrDie(task, error)
            log.trace "Task fault (3): $fault"
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
        log.trace "Finalize process > ${safeTaskName(task)}"

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
        result << "  closed: $closed\n"

        return result.toString()
    }

    class TaskProcessorInterceptor extends DataflowEventAdapter {

        final DataflowReadChannel source

        final DataflowQueue control

        final boolean singleton

        TaskProcessorInterceptor(DataflowReadChannel source, DataflowQueue control, boolean singleton) {
            this.source = source
            this.control = control
            this.singleton = singleton
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
            result[1] = messages.first()
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
            if( singleton ) {
                // -- kill the process
                control.bind(PoisonPill.instance)
            }
            else {
                // -- send a control message for each new source item to keep the process running
                control.bind(Boolean.TRUE)
            }

            return message
        }

        @Override
        Object controlMessageArrived(final DataflowProcessor processor, final DataflowReadChannel<Object> channel, final int index, final Object message) {
            // apparently auto if-guard instrumented by @Slf4j is not honoured in inner classes - add it explicitly
            if( log.isTraceEnabled() ) {
                def taskName = currentTask.get()?.name ?: name
                log.trace "<${taskName}> Control message arrived => ${message}"
            }

            super.controlMessageArrived(processor, channel, index, message)

            if( message == PoisonPill.instance ) {
                // apparently auto if-guard instrumented by @Slf4j is not honoured in inner classes - add it explicitly
                if( log.isTraceEnabled() )
                    log.trace "<${name}> Poison pill arrived; port: $index"
                closed.set(true)
                state.update { StateObj it -> it.poison() }
            }

            return message
        }

        @Override
        void afterStop(final DataflowProcessor processor) {
            // apparently auto if-guard instrumented by @Slf4j is not honoured in inner classes - add it explicitly
            if( log.isTraceEnabled() )
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
        boolean onException(final DataflowProcessor processor, final Throwable error) {
            // return `true` to terminate the dataflow processor
            handleException( error, currentTask.get() )
        }
    }
}
