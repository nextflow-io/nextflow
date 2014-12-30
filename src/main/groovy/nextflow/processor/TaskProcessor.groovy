/*
 * Copyright (c) 2013-2014, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2014, Paolo Di Tommaso and the respective authors.
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
import java.nio.file.Path
import java.nio.file.Paths
import java.util.concurrent.atomic.AtomicBoolean
import java.util.concurrent.atomic.AtomicInteger

import com.google.common.hash.HashCode
import groovy.transform.Memoized
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import groovyx.gpars.agent.Agent
import groovyx.gpars.dataflow.Dataflow
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.operator.DataflowProcessor
import groovyx.gpars.dataflow.operator.PoisonPill
import groovyx.gpars.dataflow.stream.DataflowStreamWriteAdapter
import groovyx.gpars.group.PGroup
import nextflow.Nextflow
import nextflow.Session
import nextflow.exception.MissingFileException
import nextflow.exception.MissingValueException
import nextflow.exception.ProcessException
import nextflow.exception.ProcessFailedException
import nextflow.exception.ProcessScriptException
import nextflow.executor.Executor
import nextflow.file.FileHelper
import nextflow.file.FileHolder
import nextflow.script.BaseScript
import nextflow.script.BasicMode
import nextflow.script.EachInParam
import nextflow.script.EnvInParam
import nextflow.script.FileInParam
import nextflow.script.FileOutParam
import nextflow.script.FileSharedParam
import nextflow.script.InParam
import nextflow.script.OutParam
import nextflow.script.ScriptType
import nextflow.script.SetInParam
import nextflow.script.SetOutParam
import nextflow.script.SharedParam
import nextflow.script.StdInParam
import nextflow.script.StdOutParam
import nextflow.script.TaskBody
import nextflow.script.ValueInParam
import nextflow.script.ValueOutParam
import nextflow.script.ValueSharedParam
import nextflow.util.ArrayBag
import nextflow.util.BlankSeparatedList
import nextflow.util.CacheHelper
import nextflow.util.CollectionHelper
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
abstract class TaskProcessor {

    static enum RunType {
        SUBMIT('Submitted process'),
        RETRY('Re-submitted process')

        String message;

        RunType(String str) { message=str };
    }

    static final protected String TASK_CONFIG = 'task'

    /**
     * Global count of all task instances
     */
    static final protected allCount = new AtomicInteger()

    /**
     * Unique task index number (run)
     */
    final protected indexCount = new AtomicInteger()

    /**
     * The current workflow execution session
     */
    protected final Session session

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
     */
    protected DataflowProcessor processor

    /**
     * The underlying executor which will run the task
     */
    protected final Executor executor

    /**
     * The corresponding task configuration properties, it holds the inputs/outputs
     * definition as well as other execution meta-declaration
     */
    protected final ProcessConfig config


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
     * Used to show the override warning message only the very first time
     */
    private final overrideWarnShown = new AtomicBoolean()

    /**
     * Flag set {@code true} when the processor termination has been invoked
     *
     * See {@code #checkProcessTermination}
     */
    protected volatile boolean terminated

    /**
     * Whenever the process execution is required to be blocking in order to handle
     * shared object in a thread safe manner
     */
    protected boolean blocking

    /**
     * Holds the values shared by multiple task instances
     */
    Map<SharedParam,Object> sharedObjs = [:]

    /**
     * The state is maintained by using an agent
     */
    protected Agent<StateObj> state

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
    ProcessConfig getTaskConfig() { config }

    /**
     * @return The current {@code Session} instance
     */
    Session getSession() { session }

    /**
     * @return The processor name
     */
    String getName() { name }

    /**
     * @return The {@code BaseScript} object which represents pipeline script
     */
    BaseScript getOwnerScript() { ownerScript }

    /**
     *  Define the type of script hold by the {@code #code} property
     */
    protected ScriptType getType() { taskBody.type }

    /**
     * The source fragment as entered by the user for debugging purpose
     */
    protected String getSource() { taskBody.source }

    protected Closure getCode() { taskBody.closure }

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

        if ( !code )
            throw new IllegalArgumentException("Missing 'script' attribute")


        /*
         * Normalize the input channels:
         * - at least one input channel have to be provided,
         *   if missing create an dummy 'input' set to true
         */
        log.trace "TaskConfig: ${config}"
        if( config.getInputs().size() == 0 ) {
            config.fakeInput()
        }

        final boolean hasEachParams = config.getInputs().any { it instanceof EachInParam }
        final boolean allScalarValues = config.getInputs().allScalarInputs() && !hasEachParams

        /*
         * Normalize the output
         * - even though the output may be empty, let return the stdout as output by default
         */
        if ( config.getOutputs().size() == 0 ) {
            def dummy =  allScalarValues ? Nextflow.variable() : Nextflow.channel()
            config.fakeOutput(dummy)
        }

        // the state agent
        state = new Agent<>(new StateObj(allScalarValues,name))
        state.addListener { StateObj old, StateObj obj ->
            log.trace "<$name> Process state changed to: $obj"
            if( !terminated && obj.isFinished() ) {
                terminateProcess()
                terminated = true
            }
        }

        // create the underlying dataflow operator
        createOperator()

        // register the processor
        session.taskRegister(this)

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
    protected abstract void createOperator()

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
     * Wraps the target method by a closure declaring as many arguments as many are the user declared inputs
     * object.
     *
     * @param n The number of inputs
     * @param method The method which will execute the task
     * @return The operator closure adapter
     */
    final protected Closure createCallbackWrapper( int n, Closure method ) {

        final args = []
        n.times { args << "__p$it" }

        final str = " { ${args.join(',')} -> callback([ ${args.join(',')} ]) }"
        final binding = new Binding( ['callback': method] )
        final result = (Closure)new GroovyShell(binding).evaluate (str)

        return result
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
        log.trace "Creating a new process > $name"

        def id = allCount.incrementAndGet()
        def index = indexCount.incrementAndGet()
        def task = new TaskRun(
                id: id,
                index: index,
                processor: this,
                type: type,
                config: config.createTaskConfig(),
                context: new TaskContext(this)
        )

        // setup config
        task.config.process = task.processor.name
        task.config.executor = task.processor.executor.name

        /*
         * initialize the inputs/outputs for this task instance
         */
        config.getInputs().each { InParam param ->
            if( param instanceof SetInParam )
                param.inner.each { task.setInput(it)  }
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
    final protected boolean checkCachedOrLaunchTask( TaskRun task, HashCode hash, boolean shouldTryCache, RunType runType ) {

        int tries = 0
        while( true ) {
            if( tries++ ) {
                hash = CacheHelper.defaultHasher().newHasher().putBytes(hash.asBytes()).putInt(tries).hash()
            }

            final folder = FileHelper.getWorkFolder(session.workDir, hash)
            final exists = folder.exists()
            log.trace "[${task.name}] Cacheable folder: $folder -- exists: $exists; try: $tries"

            def cached = shouldTryCache && exists && checkCachedOutput(task, folder, hash)
            if( cached )
                return false

            if( exists )
                continue

            if( !folder.mkdirs() ) {
                throw new IOException("Unable to create folder: $folder -- check file system permission")
            }

            // submit task for execution
            submitTask( task, runType, hash, folder )
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
            log.warn "[$task.name] Store dir can be used when using 'file' outputs"
            return false
        }

        if( !task.config.getStoreDir().exists() ) {
            log.trace "[$task.name] Store dir does not exists > ${task.config.storeDir} -- return false"
            // no folder -> no cached result
            return false
        }


        try {
            // -- check if all output resources are available
            collectOutputs(task)
            log.info "[skipping] Stored process > ${task.name}"

            // set the exit code in to the task object
            task.exitStatus = task.config.getValidExitStatus()[0]

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
            if( exitFile.empty() ) {
                log.trace "[$task.name] Exit file is empty > $exitFile -- return false"
                return false
            }

            def str = exitFile.text.trim()
            exitCode = str.isInteger() ? str.toInteger() : null
            if( !task.isSuccess(exitCode) ) {
                log.trace "[$task.name] Exit code is not valid > $str -- return false"
                return false
            }
        }

        /*
         * verify cached context map
         */
        TaskContext ctx = null
        def ctxFile = folder.resolve(TaskRun.CMD_CONTEXT)
        if( task.hasCacheableValues() ) {
            if( !ctxFile.exists() ) {
                log.trace "[$task.name] Context map file does not exist: $ctxFile -- return false"
                return false
            }

            ctx = TaskContext.read(this, ctxFile)
            populateSharedCtx(task, ctx)
        }

        /*
         * verify stdout file
         */
        def stdoutFile = folder.resolve( TaskRun.CMD_OUTFILE )

        try {
            // -- check if all output resources are available
            collectOutputs(task, folder, stdoutFile, ctx)

            // set the exit code in to the task object
            task.hash = hash
            task.workDir = folder
            task.stdout = stdoutFile
            if( exitCode != null ) {
                task.exitStatus = exitCode
            }
            if( ctx != null ) {
                task.context = ctx
                task.config.setContext(ctx)
                task.code?.delegate = ctx
            }

            log.info "[${task.hashLog}] Cached process > ${task.name}"

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
     * Populate the 'shared' status context map with entries
     *
     * @param task
     * @param ctx
     * @return
     */
    final protected populateSharedCtx( TaskRun task, Map ctx ) {
        task.inputs.keySet().findAll { it.class instanceof SharedParam } .each { SharedParam param ->
            sharedObjs[param.name] = ctx[param.name]
        }
    }

    /**
     * Given the script closure, that hold the code entered by the user, returns
     * the final script string to be executed
     *
     * @param code
     * @return
     */
    final protected String getScriptlet( Closure<String> code ) {
        try {
            return code.call()?.toString()
        }
        catch( Throwable e ) {
            throw new ProcessScriptException("Process script contains error(s)", e)
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
        resumeOrDie( task, error ) == ErrorStrategy.TERMINATE
    }

    /**
     *
     * @param task The {@code TaskRun} instance that raised an error
     * @param error The error object
     * @return The {@code ErrorStrategy} applied
     */
    final synchronized protected ErrorStrategy resumeOrDie( TaskRun task, Throwable error ) {
        log.trace "Handling unexpected condition for\n  task: $task\n  error [${error?.class?.name}]: ${error?.getMessage()?:error}"

        try {
            // do not recoverable error, just trow it again
            if( error instanceof Error ) throw error

            final taskErrCount = task ? task.failCount++ : 0
            final procErrCount = errorCount++
            final taskStrategy = task?.config?.getErrorStrategy()

            // when is a task level error and the user has chosen to ignore error, just report and error message
            // return 'false' to DO NOT stop the execution
            if( task && error instanceof ProcessException ) {

                // IGNORE strategy -- just continue
                if( taskStrategy == ErrorStrategy.IGNORE ) {
                    log.warn "Error running process > ${error.getMessage()} -- error is ignored"
                    task.failed = true
                    return ErrorStrategy.IGNORE
                }

                // RETRY strategy -- check that process do not exceed 'maxError' and the task do not exceed 'maxRetries'
                if( taskStrategy == ErrorStrategy.RETRY && procErrCount < task.config.getMaxErrors() && taskErrCount < task.config.getMaxRetries() ) {
                    session.execService.submit({ checkCachedOrLaunchTask( task, task.hash, false, RunType.RETRY ) } as Runnable)
                    return ErrorStrategy.RETRY
                }
            }

            // mark the task as failed
            if( task )
                task.failed = true

            // MAKE sure the error is showed only the very first time across all processes
            if( errorShown.getAndSet(true) ) {
                return ErrorStrategy.TERMINATE
            }

            def message = []
            message << "Error executing process > '${task?.name ?: name}'"
            if( error instanceof ProcessException ) {
                formatTaskError( message, error, task )
            }
            else {
                message << formatErrorCause(error)
            }
            log.error message.join('\n'), error
        }
        catch( Throwable e ) {
            // no recoverable error
            log.error("Execution aborted due to an unexpected error", e )
        }

        session.abort(error)
        return ErrorStrategy.TERMINATE
    }

    final protected formatTaskError( List message, Throwable error, TaskRun task ) {

        // compose a readable error message
        message << formatErrorCause( error )

        /*
         * task executing scriptlets
         */
        if( task?.script ) {
            // - print the executed command
            message << "Command executed:\n"
            task.script?.stripIndent()?.trim()?.eachLine {
                message << "  ${it}"
            }

            // - the exit status
            message << "\nCommand exit status:\n  ${task.exitStatus != Integer.MAX_VALUE ? task.exitStatus : '-'}"

            // - the tail of the process stdout
            message << "\nCommand output:"
            final max = 50
            def lines = task.dumpStdout(max)
            if( lines.size() == max ) {
                message << "  (more omitted..)"
            }
            else if( lines.size() == 0 ) {
                message << "  (empty)"
            }
            lines.each {
                message << "  ${task.workDir ? it.replace(task.workDir.toString()+'/','') : it }"
            }

        }
        else {
            if( source )  {
                message << "\nSource block:"
                source.stripIndent().eachLine {
                    message << "  $it"
                }
            }

        }

        if( task?.workDir )
            message << "\nWork dir:\n  ${task.workDir.toString()}"

        message << "\nTip: ${getRndTip()}"

        return message
    }

    static List tips = [
            'when you have fixed the problem you can continue the execution appending to the nextflow command line the option \'-resume\'',
            "you can try to figure out what's wrong by changing to the process work dir and showing the script file named: '${TaskRun.CMD_SCRIPT}'",
            "view the complete command output by changing to the process work dir and entering the command: 'cat ${TaskRun.CMD_OUTFILE}'",
            "you can replicate the issue by changing to the process work dir and entering the command: 'bash ${TaskRun.CMD_RUN}'"
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

        config.getOutputs().each { param ->
            def channel = param.outChannel

            if( channel instanceof DataflowQueue ) {
                log.trace "<$name> Sending a poisong pill for $param"
                channel.bind( PoisonPill.instance )
            }
            else if( channel instanceof DataflowStreamWriteAdapter ) {
                log.trace "<$name> Sending a poisong pill for $param"
                channel.bind( PoisonPill.instance )
            }
            else {
                log.trace "<$name> Poison pill is not sent over $param channel"
            }

        }
    }


    private String formatErrorCause( Throwable error ) {

        def result = new StringBuilder()
        result << '\nCaused by:\n'

        def message = error.cause?.getMessage() ?: ( error.getMessage() ?: error.toString() )
        result.append('  ').append(message).append('\n')

        result.toString()
    }



    /**
     * Bind the expected output files to the corresponding output channels
     * @param processor
     */
    synchronized protected void bindOutputs( TaskRun task ) {

        // -- creates the map of all tuple values to bind
        Map<Short,List> tuples = [:]
        config.getOutputs().each { OutParam p -> tuples.put(p.index,[]) }

        // -- collects the values to bind
        task.outputs.each { OutParam param, value ->

            switch( param ) {
            case StdOutParam:
                log.trace "Process $name > normalize stdout param: $param"
                value = value instanceof Path ? value.text : value?.toString()

            case FileOutParam:
            case ValueOutParam:
                log.trace "Process $name > collecting out param: ${param} = $value"
                tuples[param.index].add( value )
                break

            default:
                throw new IllegalArgumentException("Illegal output parameter type: $param")
            }
        }

        // -- bind out the collected values
        config.getOutputs().each { param ->
            def list = tuples[param.index]
            if( list == null ) throw new IllegalStateException()

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

    protected void bindOutParam( OutParam param, def values ) {
        log.trace "<$name> Binding param $param with $values"
        def x = values.size() == 1 ? values[0] : values
        processor.bindOutput( param.index, x )
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

        /*
         * Shared objects behave as accumulators
         * Copying back updated values from context map to buffer map so that can be accessed in the next iteration
         */
        log.trace "<${name}> collecting sharedObjs: $sharedObjs"
        sharedObjs?.keySet()?.each { param ->

            switch(param) {
                case ValueSharedParam:
                    sharedObjs[param] = context.get(param.name)
                    break

                case FileSharedParam:
                    // the 'sharedObjs' is updated in the 'beforeRun' method
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
     * @param param The declared {@code StdOutParam} object
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

        def all = []
        def fileParam = param as FileOutParam
        // type file parameter can contain a multiple files pattern separating them with a special character
        def entries = param.getFilePatterns(context)

        // for each of them collect the produced files
        entries.each { String pattern ->
            List<Path> result=null
            if( FileHelper.isGlobPattern(pattern) ) {
                result = executor.collectResultFile(workDir, pattern, task.name, param)
                // filter the inputs
                if( !fileParam.includeInputs ) {
                    result = filterByRemovingStagedInputs(task, result)
                    log.trace "Process ${task.name} > after removing staged inputs: ${result}"
                }
            }
            else {
                def file = workDir.resolve(pattern)
                if( file.exists() )
                    result = [file]
            }

            if( !result )
                throw new MissingFileException("Missing output file(s): '$pattern' expected by process: ${task.name}")

            all.addAll(result)
        }

        task.setOutput( param, all.size()==1 ? all[0] : all )

    }

    protected void collectOutValues( TaskRun task, ValueOutParam param, Map ctx ) {

        // look into the task inputs value for an *ValueInParam* entry
        // having the same *name* as the requested output name
        if( !ctx.containsKey(param.name) ) {
            throw new MissingValueException("Missing value declared as output paramter: ${param.name}")
        }

        // bind the value
        def val = ctx.get(param.name)
        task.setOutput( param, val )
        log.trace "Collecting param: ${param.name}; value: ${val}"

    }


    /**
     * Given a list of {@code Path} removes all the hidden file i.e. the ones which names starts with a dot char
     * @param files A list of {@code Path}
     * @return The result list not containing hidden file entries
     */
    protected List<Path> filterByRemovingHiddenFiles( List<Path> files ) {
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
    protected List<Path> filterByRemovingStagedInputs( TaskRun task, List<Path> files ) {

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
                result['PATH'] =  "${session.binDir}:${result['PATH']}".toString()
            }
            else {
                result['PATH'] = "${session.binDir}:\$PATH".toString()
            }
        }

        return result
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

        if( input instanceof Path ) {
            return new FileHolder(input)
        }

        def result = Nextflow.tempFile(altName)
        def source = input?.toString() ?: ''
        result.text = source
        return new FileHolder(source, result)
    }


    protected List<FileHolder> normalizeInputToFiles( Object obj, int count ) {

        Collection allItems = obj instanceof Collection ? obj : [obj]
        def len = allItems.size()

        // use a bag so that cache hash key is not affected by file entries order
        def files = new ArrayBag(len)
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
    protected List<FileHolder> expandWildcards( String name, List<FileHolder> files ) {
        assert files != null

        // use an unordered so that cache hash key is not affected by file entries order
        final result = new ArrayBag()
        if( files.size()==0 ) { return result }

        if( !name || name == '*' ) {
            return files
        }

        // no wildcards in the file name
        else if( !name.contains('*') && !name.contains('?') ) {
            /*
             * The name do not contain any wildcards *BUT* when multiple files are provide
             * it is managed like having a 'start' at the end of the file name
             */
            if( files.size()>1 ) {
                name += '*'
            }
            else {
                // just return that name
                result << files[0].withName(name)
                return result
            }
        }

        /*
         * The star wildcard: when a single item is provided, it is simply ignored
         * When a collection of files is provided, the name is expanded to the index number
         */
        if( name.contains('*') ) {
            if( files.size()>1 ) {
                def count = 1
                files.each { FileHolder holder ->
                    def stageName = name.replace('*', (count++).toString())
                    result << holder.withName(stageName)
                }
            }
            else {
                // there's just one value, remove the 'star' wildcards
                def stageName = name.replace('*','')
                result << files[0].withName(stageName)
            }
        }

        /*
         * The question mark wildcards *always* expand to an index number
         * as long as are the number of question mark characters
         */
        else if( name.contains('?') ) {
            def count = 1
            files.each { FileHolder holder ->
                String match = (name =~ /\?+/)[0]
                def replace = (count++).toString().padLeft(match.size(), '0')
                def stageName = name.replace(match, replace)
                result << holder.withName(stageName)
            }

        }

        // not a valid condition
        else {
            throw new IllegalStateException("Invalid file expansion for name: '$name'")
        }

        return result
    }


    /**
     * Given a map holding variables key-value pairs, create a script fragment
     * exporting the required environment variables
     */
    static String bashEnvironmentScript( Map<String,String> environment ) {

        final script = []
        environment.each { name, value ->
            if( name ==~ /[a-zA-Z_]+[a-zA-Z0-9_]*/ ) {
                script << "export $name=\"$value\""
            }
            else {
                log.trace "Illegal environment variable name: '${name}'"
            }
        }
        script << ''

        return script.join('\n')
    }


    protected decodeInputValue( InParam param, List values ) {

        def val = values[ param.index ]
        if( param.mapIndex != -1 ) {
            def list
            if( val instanceof Map )
                list = val.values()

            else if( val instanceof Collection )
                list = val
            else
                list = [val]

            try {
                return list[param.mapIndex]
            }
            catch( IndexOutOfBoundsException e ) {
                throw new ProcessException(e)
            }
        }

        return val
    }

    /**
     * Create the {@code TaskDef} data structure and initialize the task execution context
     * with the received input values
     *
     * @param values
     * @return
     */
    final protected TaskRun setupTask(List values) {
        log.trace "Setup new process > $name"

        // -- map the inputs to a map and use to delegate closure values interpolation
        final secondPass = [:]
        final task = createTaskRun()

        int count = makeTaskContextStage1(task, secondPass, values)
        makeTaskContextStage2(task, secondPass, count)

        return task
    }

    final protected int makeTaskContextStage1( TaskRun task, Map secondPass, List values ) {

        final contextMap = task.context
        final firstRun = task.index == 1
        int count = 0

        task.inputs.keySet().each { InParam param ->

            // add the value to the task instance
            def val = decodeInputValue(param,values)

            switch(param) {
                case EachInParam:
                case ValueInParam:
                    contextMap.put( param.name, val )
                    break

                case FileInParam:
                    secondPass[param] = val
                    return // <-- leave it, because we do not want to add this 'val' at this stage

                case FileSharedParam:
                    def fileParam = param as FileSharedParam
                    if( firstRun ) {
                        def normalized = normalizeInputToFiles(val,count)
                        if( normalized.size() > 1 )
                            throw new IllegalStateException("Cannot share multiple files")

                        val = expandWildcards( fileParam.filePattern, normalized )
                        count += val.size()
                        // track this obj
                        sharedObjs[(SharedParam)param] = val
                    }
                    else {
                        val = sharedObjs[(SharedParam)param]
                    }

                    contextMap.put( fileParam.name, singleItemOrList(val) )
                    break

                case ValueSharedParam:
                    if( firstRun )
                        sharedObjs[(SharedParam)param] = val
                    else
                        val = sharedObjs[(SharedParam)param]

                    contextMap.put( param.name, val )
                    break

                case StdInParam:
                case EnvInParam:
                    // nothing to do
                    break

                default:
                    log.debug "Unsupported input param type: ${param?.class?.simpleName}"
            }

            // add the value to the task instance context
            task.setInput(param, val)
        }

        return count
    }

    final protected void makeTaskContextStage2( TaskRun task, Map secondPass, int count ) {

        final ctx = task.context

        // -- all file parameters are processed in a second pass
        //    so that we can use resolve the variables that eventually are in the file name
        secondPass.each { FileInParam param, val ->

            def fileParam = param as FileInParam
            def normalized = normalizeInputToFiles(val,count)
            def resolved = expandWildcards( fileParam.getFilePattern(ctx), normalized )
            ctx.put( param.name, singleItemOrList(resolved) )
            count += resolved.size()

            // add the value to the task instance context
            task.setInput(param, resolved)
        }

        // -- set the delegate map as context ih the task config
        //    so that lazy directives will be resolved against it
        task.config.setContext(ctx)

        // -- set the `task.config` object in the context
        if( !ctx.containsKey(TASK_CONFIG) ) {
            ctx.put(TASK_CONFIG, task.config)
        }
        else if( !overrideWarnShown.getAndSet(true) ) {
            log.warn "Process $name overrides reserved variable `task`"
        }

        // -- initialize the task code to be executed
        task.code = this.code.clone() as Closure
        task.code.delegate = task.context
        task.code.setResolveStrategy(Closure.DELEGATE_ONLY)

        // Important!
        // when the task is implemented by a script string
        // Invokes the closure which return the script whit all the variables replaced with the actual values
        if( type == ScriptType.SCRIPTLET ) {
            task.script = getScriptlet(task.code)
        }
    }

    final protected void makeTaskContextStage3( TaskRun task, HashCode hash, Path folder ) {

        // set hash-code & working directory
        task.hash = hash
        task.workDir = folder
        task.config.workDir = folder

    }

    /**
     * Execute the specified task shell script
     *
     * @param script The script string to be execute, e.g. a BASH script
     * @return {@code TaskDef}
     */
    final protected void submitTask( TaskRun task, RunType runType, HashCode hash, Path folder ) {
        log.trace "[${task.name}] actual run folder: ${task.workDir}"

        makeTaskContextStage3(task, hash, folder)

        // add the task to the collection of running tasks
        session.dispatcher.submit(task, blocking, runType.message)
    }


    /**
     * Finalize the task execution, checking the exit status
     * and binding output values accordingly
     *
     * @param task The {@code TaskRun} instance to finalize
     */
    @PackageScope
    final void finalizeTask( TaskRun task ) {
        log.trace "finalizing process > ${task.name} -- $task"

        def strategy = null
        try {
            if( task.type == ScriptType.SCRIPTLET ) {
                if( task.exitStatus == Integer.MAX_VALUE )
                    throw new ProcessFailedException("Process '${task.name}' terminated for an unknown reason -- Likely it has been terminated by the external system")

                if ( !task.isSuccess() )
                    throw new ProcessFailedException("Process '${task.name}' terminated with an error exit status")
            }

            // -- verify task exist status
            if( task.error )
                throw new ProcessFailedException("Process '${task.name}' failed", task.error)

            // -- if it's OK collect results and finalize
            collectOutputs(task)

            // save the context map for caching purpose
            // only the 'cache' is active and
            if( isCacheable() && task.hasCacheableValues() && task.context != null ) {
                def target = task.workDir.resolve(TaskRun.CMD_CONTEXT)
                if( task.context.get(TASK_CONFIG) instanceof TaskConfig )
                    task.context.remove(TASK_CONFIG)
                task.context.save(target)
            }

        }
        catch ( Throwable error ) {
            strategy = resumeOrDie(task, error)
        }

        // -- finalize the task
        if( strategy != ErrorStrategy.RETRY )
            finalizeTask0(task)
    }

    /**
     * Whenever the process can be cached
     */
    protected boolean isCacheable() {
        session.cacheable && config.cacheable
    }

    protected boolean isResumable() {
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
        }

        // increment the number of processes executed
        state.update { StateObj it -> it.incCompleted() }
    }

    protected void terminateProcess() {
        log.debug "<${name}> Sending poison pills and terminating process"
        sendPoisonPill()
        session.taskDeregister(this)
        processor.terminate()
    }

}

