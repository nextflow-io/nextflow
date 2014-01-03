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

import java.nio.file.Files
import java.nio.file.Path
import java.util.concurrent.atomic.AtomicInteger
import java.util.concurrent.locks.ReentrantLock

import embed.com.google.common.hash.HashCode
import groovy.transform.PackageScope
import groovy.transform.Synchronized
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.Dataflow
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowVariable
import groovyx.gpars.dataflow.operator.DataflowProcessor
import groovyx.gpars.dataflow.operator.PoisonPill
import groovyx.gpars.dataflow.stream.DataflowStreamWriteAdapter
import groovyx.gpars.group.PGroup
import nextflow.Nextflow
import nextflow.Session
import nextflow.ast.AstNodeToScriptVisitor
import nextflow.exception.ProcessFailedException
import nextflow.exception.MissingFileException
import nextflow.exception.MissingValueException
import nextflow.exception.ProcessException
import nextflow.exception.ProcessScriptException
import nextflow.executor.AbstractExecutor
import nextflow.script.BaseScript
import nextflow.script.ScriptType
import nextflow.util.BlankSeparatedList
import nextflow.util.CacheHelper
import nextflow.util.FileHelper
import org.apache.commons.lang.SerializationUtils

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
abstract class TaskProcessor {

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
    protected final BaseScript ownerScript

    /**
     * Gpars thread pool
     */
    protected PGroup group = Dataflow.retrieveCurrentDFPGroup()

    /**
     * The processor descriptive name
     */
    protected String name

    /**
     * The closure wrapping the script to be executed
     */
    protected Closure code

    /**
     *  Define the type of script hold by the {@code #code} property
     */
    protected ScriptType type = ScriptType.SCRIPTLET

    /**
     * The source fragment as entered by the user for debugging purpose
     */
    protected String source

    /**
     * The corresponding {@code DataflowProcessor} which will receive and
     * manage accordingly the task inputs
     */
    protected DataflowProcessor processor

    /**
     * Whenever all inputs for this task are *scalar* value, i.e. simple data types and not collections of values
     */
    protected volatile boolean allScalarValues

    /**
     * The underlying executor which will run the task
     */
    protected final AbstractExecutor executor

    /**
     * The corresponding task configuration properties, it holds the inputs/outputs
     * definition as well as other execution meta-declaration
     */
    protected final TaskConfig taskConfig

    /**
     * Lock to protected against race-condition the folder creation phase
     */
    private static final folderLock = new ReentrantLock(true)

    /**
     * Used for framework generated task names
     */
    @Deprecated
    private static final AtomicInteger tasksCount = new AtomicInteger()

    /*
     * Internally used, random number generator
     */
    private final random = new Random()

    /**
     * Count the number errors showed
     */
    private final errorCount = new AtomicInteger()

    /**
     * Count how many times the process finalization method has been invoked
     *
     * See {@code #finalizeTask0()}
     * See {@code #checkProcessTermination}
     */
    protected final finalizeCount = new AtomicInteger()

    /**
     * Count how many times this process has been launched
     */
    protected final instanceCount = new AtomicInteger()

    /**
     * Flat set {@code true} when the processor receive a poison pill message (to stop it)
     */
    protected volatile boolean receivedPoisonPill

    /**
     * Flag set {@code true} when the processor termination has been invoked
     *
     * See {@code #checkProcessTermination}
     */
    protected volatile boolean terminated

    /**
     * Whener the process execution is required to be blocking in order to handle
     * shared object in a thread safe manner
     */
    protected boolean blocking

    /**
     * Holds the values shared by multiple task instances
     */
    Map<SharedParam,Object> sharedObjs = [:]

    /* for testing purpose - do not remove */
    protected TaskProcessor() { }

    /**
     * Create and initialize the processor object
     *
     * @param executor
     * @param session
     * @param script
     * @param taskConfig
     * @param taskBlock
     */
    TaskProcessor( AbstractExecutor executor, Session session, BaseScript script, TaskConfig taskConfig, Closure taskBlock ) {
        assert executor
        assert session
        assert script
        assert taskBlock

        this.executor = executor
        this.session = session
        this.ownerScript = script
        this.taskConfig = taskConfig
        this.code = taskBlock

        /*
         * set the task name
         */
        if( taskConfig.name ) {
            this.name = taskConfig.name
        }
        else {
            // generate the processor name if not specified
            this.name = "task_${tasksCount.incrementAndGet()}"
            taskConfig.name = this.name
        }
    }

    /**
     * @return The {@code TaskConfig} object holding the task configuration properties
     */
    TaskConfig getTaskConfig() { taskConfig }

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

        if ( !code ) {
            throw new IllegalArgumentException("Missing 'script' attribute")
        }


        /*
         * Normalize the input channels:
         * - at least one input channel have to be provided,
         *   if missing create an dummy 'input' set to true
         */
        log.trace "TaskConfig: ${taskConfig}"
        if( taskConfig.inputs.size() == 0 ) {
            taskConfig.noInput()
        }

        allScalarValues = !taskConfig.inputs.channels.any { !(it instanceof DataflowVariable) }

        /*
         * Normalize the output
         * - even though the output may be empty, let return the stdout as output by default
         */
        if ( taskConfig.outputs.size() == 0 ) {
            def dummy =  allScalarValues ? Nextflow.val() : Nextflow.channel()
            taskConfig.stdout(dummy)
        }

        // create the underlying dataflow operator
        createOperator()

        // register the processor
        session.taskRegister()

        /*
         * When there is a single output channel, return let returns that item
         * otherwise return the list
         */
        def result = taskConfig.outputs.channels
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
    def getShebangLine() {
        assert taskConfig.shell, "Missing 'shell' property for process: $name"

        def shell = taskConfig.shell
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
     * Remove extra leading and trailing whitespace and newlines chars,
     * also if the script does not start with a {@code shebang} line,
     * add the default by using the current {@code #shell} attribute
     */
    def String normalizeScript(String script) {
        assert script != null

        def result = new StringBuilder()
        result << script.stripIndent().trim()
        result << '\n'

        if( result[0] != '#' || result[1] != '!') {
            result.insert(0, shebangLine +'\n')
        }

        return result.toString()
    }

    /**
     * Given the task script extract the top *she-bang* interpreter declaration removing the {@code #!} characters.
     * @param script The script to be executed
     *
     * @return The interpreter as defined in the she-bang declaration, for example {@code /usr/bin/env perl}
     */
    def String fetchInterpreter( String script ) {
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
     * @return The new newly created {@code TaskRun{
     */
    final protected TaskRun createTaskRun() {
        log.trace "Creating a new process > $name"

        def id = allCount.incrementAndGet()
        def index = indexCount.incrementAndGet()
        new TaskRun(id: id, index: index, name: "$name ($index)", processor: this, type: type )
    }

    final protected getHashLog( HashCode hash ) {
        def str = hash.toString()
        def result = new StringBuilder()
        result << str[0]
        result << str[1]
        result << '/'
        for( int i=2; i<8 && i<str.size(); i++ ) {
            result << str[i]
        }
        return result.toString()
    }

    /**
     * Create a unique-folder where the task will be executed
     *
     * @param folder
     * @param hash
     * @return
     */
    final protected Path createTaskFolder( Path folder, HashCode hash  ) {

        folderLock.lock()
        try {

            if( folder.exists() ) {

                // find another folder name that does NOT exist
                while( true ) {
                    hash = CacheHelper.hasher( [hash.asInt(), random.nextInt() ] ).hash()
                    folder = FileHelper.getWorkFolder(session.workDir, hash)
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

    /**
     * Check whenever the outputs for the specified task already exist
     *
     * @param task The task instance
     * @param folder The folder where the outputs are stored (eventually)
     * @return {@code true} when all outputs are available, {@code false} otherwise
     */
    final boolean checkCachedOutput(TaskRun task, Path folder, HashCode hash) {
        if( !folder.exists() ) {
            log.trace "[$task.name] Cached folder does not exists > $folder -- return false"
            // no folder -> no cached result
            return false
        }

        // check if exists the task exit code file
        def exitCode = null
        def exitFile = folder.resolve(TaskRun.CMD_EXIT)
        if( task.type == ScriptType.SCRIPTLET ) {
            if( exitFile.empty() ) {
                log.trace "[$task.name] Exit file is empty > $exitFile -- return false"
                return false
            }

            def exitValue = exitFile.text.trim()
            exitCode = exitValue.isInteger() ? exitValue.toInteger() : null
            if( exitCode == null || !(exitCode in taskConfig.validExitStatus) ) {
                log.trace "[$task.name] Exit code is not valid > $exitValue -- return false"
                return false
            }
        }

        /*
         * verify cached context map
         */
        def ctxMap = null
        def ctxFile = folder.resolve(TaskRun.CMD_CONTEXT)
        def outCount = task.getOutputsByType(ValueOutParam).size()
        if( outCount ) {
            if( !ctxFile.exists() ) {
                log.trace "[$task.name] Contexy map file does not exist: $ctxFile -- return false"
                return false
            }
            ctxMap = readContextMap(ctxFile)
        }

        /*
         * verify stdout file
         */
        def stdoutFile = folder.resolve( TaskRun.CMD_OUTFILE )

        try {
            // -- check if all output resources are available
            collectOutputs(task, folder, stdoutFile, ctxMap)
            log.info "[${getHashLog(hash)}] Cached process > ${task.name}"

            // set the exit code in to the task object
            task.workDirectory = folder
            task.stdout = stdoutFile
            if( exitCode != null ) {
                task.exitStatus = exitCode
            }
            if( task.code && ctxMap ) {
                task.code.delegate = ctxMap
            }

            // -- now bind the results
            finalizeTask0(task)
            return true
        }
        catch( MissingFileException | MissingValueException e ) {
            log.debug "[$task.name] Missed cache > ${e.getMessage()} -- folder: $folder"
            task.exitStatus = Integer.MAX_VALUE
            task.workDirectory = null
            return false
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

        // when is a task level error and the user has chosen to ignore error, just report and error message
        // return 'false' to DO NOT stop the execution
        if( error instanceof ProcessException && taskConfig.errorStrategy == ErrorStrategy.IGNORE ) {
            log.warn "Error running process > ${error.getMessage()} -- error is ignored"
            return false
        }

        // make sure the error is showed only the very first time
        if( errorCount.getAndIncrement()>0 ) {
            return true
        }

        def message = []
        message << "Error executing process > '${task?.name ?: name}'"
        if( error instanceof ProcessException ) {
            formatTaskError( message, error, task )
        }
        else {
            message << formatErrorCause( error )
            message << "Tip: check the log file '.nextflow.log' for more details"
        }

        log.error message.join('\n')
        log.debug "Error details", error

        session.abort()
        return true
    }

    final protected formatTaskError( List message, Throwable error, TaskRun task ) {

        // compose a readable error message

        /*
         * task executing scriptlets
         */
        if( task?.script ) {
            // print the executed command
            message << "Command executed:\n"
            task.script.eachLine {
                message << "  $it"
            }

            message << "\nCommand exit status:\n  ${task.exitStatus != Integer.MAX_VALUE ? task.exitStatus : '-'}"

            message << "\nCommand output:"
            task.stdout?.eachLine {
                message << "  $it"
            }

        }
        else {
            if( source )  {
                message << "\nSource block:"
                source.stripIndent().eachLine {
                    message << "  $it"
                }
            }

            message << formatErrorCause( error )
        }


        if( task?.workDirectory )
            message << "\nWork dir:\n  ${task.workDirectory.toString()}"

        message << "\nTip: when you have fixed the problem you may continue the execution appending to the nextflow command line the '-resume' option"

        return message
    }


    /**
     * Send a poison pill over all the outputs channel
     */
    final protected synchronized void sendPoisonPill() {
        log.trace "Forwarding Poison-pill for process > ${name}"

        taskConfig.outputs.eachParam { name, channel ->
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

    private String getCodeScript( TaskRun task ) {

        try {
            def node = task.code.metaClass.classNode.getDeclaredMethods("doCall")[0].code
            def writer = new StringWriter()
            node.visit( new AstNodeToScriptVisitor(writer) )
            writer.toString()
        }
        catch( Throwable e ) {
            log.debug "Unable to obtain code for process: ${task.name}"
            return null
        }

    }


    private String formatErrorCause( Throwable error ) {

        def result = new StringBuilder()
        result << '\nError cause:'

        def message = error.cause?.getMessage() ?: ( error.getMessage() ?: error.toString() )
        result.append('\n  ').append(message)

        result.toString()
    }


    /**
     * Bind the expected output files to the corresponding output channels
     * @param processor
     */
    synchronized protected void bindOutputs( TaskRun task ) {

        // -- bind each produced file to its own channel
        task.outputs.eachWithIndex { OutParam param, value, index ->

            switch( param ) {
            case StdOutParam:
                log.trace "Process $name > Binding '$value' to stdout"
                processor.bindOutput(index, value instanceof Path ? value.text : value?.toString())
                break

            case FileOutParam:
                log.trace "Process $name > Binding file: '$value' to '${param.name}'"
                if( value instanceof Collection && (param as FileOutParam).flat ) {
                    value.each { processor.bindOutput(index, it) }
                }
                else {
                    processor.bindOutput(index, value)
                }
                break;

            case ValueOutParam:
                log.trace "Process $name > Binding value: '$value' to '${param.name}'"
                processor.bindOutput(index, value)
                break

            default:
                throw new IllegalArgumentException("Illegal output parameter type: ${param.class.simpleName}")
            }
        }

        // -- finally prints out the task output when 'echo' is true
        if( taskConfig.echo ) {
            task.echoStdout()
        }
    }

    final protected void collectOutputs( TaskRun task ) {
        collectOutputs( task, task.workDirectory, task.@stdout, task.code?.delegate )
    }

    /**
     * Once the task has completed this method is invoked to collected all the task results
     *
     * @param task
     */
    final protected void collectOutputs( TaskRun task, Path workDir, def stdout, Map context ) {

        task.outputs.keySet().each { OutParam param ->

            switch( param ) {
                case StdOutParam:
                    collectStdOut(task, param, stdout)
                    break

                case FileOutParam:
                    collectOutFiles(task, param, workDir)
                    break

                case ValueOutParam:
                    collectOutValues(task, param, context)
                    break

                default:
                    throw new IllegalArgumentException("Illegal output parameter: ${param.class.simpleName}")

            }
        }

        /*
         * Shared objects behave as accumulators
         * Copying back updated values from context map to buffer map so that can be accessed in the next iteration
         *
         */
        sharedObjs?.keySet() .each { param ->

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


    protected void collectOutFiles( TaskRun task, FileOutParam param, Path workDir ) {

        def all = []
        def fileParam = param as FileOutParam
        // type file parameter can contain a multiple files pattern separating them with a special character
        def entries = fileParam.separatorChar ? fileParam.name.split(/\${fileParam.separatorChar}/) : [fileParam.name]
        // for each of them collect the produced files
        entries.each { String pattern ->
            def result = executor.collectResultFile(workDir, pattern, task.name)
            log.debug "Process ${task.name} > collected outputs for pattern '$pattern': $result"

            if( result instanceof List ) {
                // filter the result collection
                if( pattern.startsWith('*') && !fileParam.includeHidden ) {
                    result = filterByRemovingHiddenFiles(result)
                    log.trace "Process ${task.name} > after removing hidden files: ${result}"
                }

                // filter the inputs
                if( !fileParam.includeInputs ) {
                    result = filterByRemovingStagedInputs(task, result)
                    log.trace "Process ${task.name} > after removing staged inputs: ${result}"
                }

                all.addAll((List) result)
            }

            else if( result ) {
                all.add(result)
            }
        }

        task.setOutput( param, all.size()==1 ? all[0] : all )

    }

    protected void collectOutValues( TaskRun task, ValueOutParam param, Map ctx ) {

        // look into the task inputs value for an *ValueInParam* entry
        // having the same *name* as the requested output name
        if( !ctx.containsKey(param.name) ) {
            throw new MissingValueException("Illegal output val parameter: ${param.name}")
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
        def List<Path> allStaged = task.getStagedInputs()
        def List<String>  allInputNames = allStaged.collect { it.getName() }

        files.findAll { !allInputNames.contains(it.getName()) }

    }

    /**
     * @return The map holding the shell environment variables for the task to be executed
     */
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
    protected FileHolder normalizeInputToFile( Object input, String altName, FileSpec fileSpec = null ) {

        if( input instanceof Path ) {
            checkSpec(input,fileSpec)
            return new FileHolder(input)
        }

        def result = ownerScript.tempFile(altName)
        def source = input?.toString() ?: ''
        result.text = source
        return new FileHolder(source, result)
    }

    protected void checkSpec(Path path, FileSpec fileSpec) {

        if( !path.exists() && fileSpec?.create ) {
            if( fileSpec.file )
                Files.createFile(path)
            else
                Files.createDirectory(path)
        }

    }

    protected List<FileHolder> normalizeInputToFiles( Object obj, int count, FileSpec fileSpec = null ) {

        def files = (obj instanceof Collection ? obj : [obj]).collect {
            normalizeInputToFile(it, "input.${++count}", fileSpec)
        }

        return files
    }

    protected singleItemOrList( List<FileHolder> items ) {
        assert items
        if( items.size() == 1 ) {
            return items[0].stagePath
        }
        else {
            return new BlankSeparatedList( items *. stagePath )
        }
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

        final result = []
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
                result << holder.withName( stageName )
            }

        }

        // not a valid condition
        else {
            throw new IllegalStateException("Invalid file expansion for name: '$name'")
        }

        return result
    }


    @Deprecated
    protected List<FileHolder> expandWildcards( String filePattern, FileHolder... files ) {
        expandWildcards( filePattern, files as List )
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



    /**
     * Execute the specified task shell script
     *
     * @param script The script string to be execute, e.g. a BASH script
     * @return {@code TaskDef}
     */
    final protected void submitTask( TaskRun task ) {
        // add the task to the collection of running tasks
        session.dispatcher.submit(task, blocking)
    }


    /**
     * Map used to delegate variable resolution to script scope
     */
    @Slf4j
    static class DelegateMap implements Map {

        @Delegate
        private Map<String,Object> holder

        private BaseScript script

        DelegateMap(BaseScript script) {
            this.script = script
            this.holder = [:]
        }

        DelegateMap(BaseScript script, Map holder) {
            assert holder != null
            this.script = script
            this.holder = holder
        }

        @Override
        public Object get(Object property) {

            def result = null
            if( holder.containsKey(property) ) {
                return holder.get(property)
            }
            // todo -- verify if it could better using "script.getBinding().getVariable()"
            else if ( script && script.hasProperty(property)) {
                return script.getProperty(property?.toString())
            }

            throw new MissingPropertyException("Unknown variable '$property' -- Make sure you didn't misspell it or define somewhere in the script before use it")

        }

        @Override
        public put(String property, Object newValue) {
            holder.put(property, newValue)
        }

        public Map getHolder() {
            holder
        }
    }


    /**
     * Finalize the task execution, checking the exit status
     * and binding output values accordingly
     *
     * @param task The {@code TaskRun} instance to finalize
     */
    @PackageScope
    final void finalizeTask( TaskRun task ) {
        log.trace "finalizing process > ${task.name}"
        try {
            // verify task exist status
            if( task.type == ScriptType.GROOVY ) {
                if( task.error ) {
                    throw new ProcessFailedException("Process '${task.name}' failed", task.error)
                }
            }

            else {
                boolean success = (task.exitStatus in taskConfig.validExitStatus)
                if ( !success ) {
                    throw new ProcessFailedException("Process '${task.name}' failed")
                }
            }

            // if it's OK collect results and finalize
            collectOutputs(task)

            // save the context map for caching purpose
            // only the 'cache' is active and
            def cacheable = session.cacheable && task.processor.taskConfig.cacheable
            if( cacheable && task.getOutputsByType(ValueOutParam).size() ) {
                saveContextMap( (DelegateMap) task.code.delegate, task.getCmdContextFile() )
            }

        }
        catch ( Throwable error ) {
            handleException(error, task)
        }
        finally {
            finalizeTask0(task)
        }
    }

    /**
     * Finalize the task execution, checking the exit status
     * and binding output values accordingly
     *
     * @param task The {@code TaskRun} instance to finalize
     * @param producedFiles The map of files to be bind the outputs
     */
    private void finalizeTask0( TaskRun task ) {
        log.debug "Finalize process > ${task.name}"

        // -- bind output (files)
        if( task.canBind ) {
            bindOutputs(task)
        }

        checkProcessTermination(true)
    }


    @Synchronized
    protected boolean checkProcessTermination( boolean isFinalize = false ) {

        if( terminated ) {
            return true
        }

        def created = instanceCount.get()
        def finalized = isFinalize ? finalizeCount.incrementAndGet() : finalizeCount.get()
        // log.debug "Finalizing task > task: ${name}; finalize: $isFinalize; allScalarValues: ${allScalarValues}; receivedPoisonPill: ${receivedPoisonPill}; instancesCount: ${tot}; finalizeCount: ${count} "

        def done = allScalarValues || ( receivedPoisonPill && created == finalized )
        if( done ) {
            log.debug "Finalizing process > ${name} -- isFinalize: $isFinalize"
            sendPoisonPill()
            session.taskDeregister()
            processor.terminate()
            terminated = true
        }

        return done
    }


    protected void saveContextMap( DelegateMap map, Path contextFile ) {
        contextFile.bytes = SerializationUtils.serialize( map.getHolder() )
    }

    protected DelegateMap readContextMap( Path contextFile ) {
        def map = (Map)SerializationUtils.deserialize( contextFile.bytes )
        new DelegateMap(ownerScript,map)
    }

}

