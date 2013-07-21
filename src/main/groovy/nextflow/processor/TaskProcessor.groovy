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
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.Dataflow
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowVariable
import groovyx.gpars.dataflow.operator.DataflowProcessor
import groovyx.gpars.dataflow.operator.PoisonPill
import groovyx.gpars.dataflow.stream.DataflowStreamWriteAdapter
import groovyx.gpars.group.PGroup
import nextflow.Session
import nextflow.exception.MissingFileException
import nextflow.exception.TaskValidationException
import nextflow.executor.AbstractExecutor
import nextflow.script.BaseScript
import nextflow.util.CacheHelper
import nextflow.util.FileHelper
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
abstract class TaskProcessor {

    protected int index

    /**
     * The current workflow execution session
     */
    protected final Session session

    /**
     * The script object which defines this task
     */
    protected final BaseScript ownerScript

    protected PGroup group = Dataflow.retrieveCurrentDFPGroup()

    /**
     * The processor descriptive name
     */
    protected String name

    /**
     * The closure wrapping the script to be executed
     */
    protected Closure code

    protected TaskRun lastRunTask

    protected DataflowProcessor processor

    protected allScalarValues

    protected final AbstractExecutor executor

    protected final TaskConfig taskConfig

    private static final creationLock = new ReentrantLock(true)

    private static final folderLock = new ReentrantLock(true)

    private static final AtomicInteger tasksCount = new AtomicInteger()

    private final random = new Random()

    private Boolean errorShown = Boolean.FALSE

    private final errorLock = new ReentrantLock(true)


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

        executor.taskConfig = taskConfig

        if( taskConfig.name ) {
            this.name = taskConfig.name
        }
        else {
            /*
             * generate the processor name if not specified
             */
            this.name = "task_${tasksCount.incrementAndGet()}"
            taskConfig.name = this.name
        }
    }


    TaskConfig getTaskConfig() {
        taskConfig
    }

    Session getSession() { session }

    String getName() { name }

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
            taskConfig.input('$':true)
        }

        allScalarValues = !taskConfig.inputs.values().any { !(it instanceof DataflowVariable) }

        /*
         * Normalize the output
         * - even though the output may be empty, let return the stdout as output by default
         */
        if ( taskConfig.outputs.size() == 0 ) {
            taskConfig.output('-')
        }

        /*
         * bind the outputs to the script scope
         */
        if( ownerScript ) {
            taskConfig.outputs.each { String name, Object channel ->
                if( name != '-' ) { ownerScript.setProperty(name, channel) }
            }
        }

        createOperator()

        /*
         * When there is a single output channel, return let returns that item
         * otherwise return the list
         */
        def result = taskConfig.outputs.values()
        return result.size() == 1 ? result[0] : result
    }

    protected abstract void createOperator()

    def getShebangLine() {

        def shell = taskConfig['shell']
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

    final protected TaskRun createTaskRun() {

        TaskRun result = null
        creationLock.lock()
        try {
            def num = session.tasks.size()
            result = new TaskRun(processor: this, id: num, status: TaskRun.Status.PENDING, index: ++index, name: "$name ($index)" )
            session.tasks.put( this, result )
        }
        finally {
            creationLock.unlock()
        }

        return result
    }


    final protected File createTaskFolder( File folder, HashCode hash  ) {

        folderLock.lock()
        try {

            if( folder.exists() ) {

                // find another folder name that does NOT exist
                while( true ) {
                    hash = CacheHelper.hasher( [hash.asInt(), random.nextInt() ] ).hash()
                    folder = FileHelper.createWorkFolder(session.workDir, hash)
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
        if( exitCode == null || !(exitCode in taskConfig.validExitCodes) ) {
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
            if( taskConfig.echo ) {
                def out = executor.getStdOutFile(task)
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
        if( error instanceof TaskValidationException && taskConfig.errorStrategy == ErrorStrategy.IGNORE ) {
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

        // -- when all input values are 'scalar' (not queue or broadcast stream)
        //    stops after the first run
        if( allScalarValues ) {
            log.debug "Finalize > ${name} terminates since all values are scalar"
            // send a poison pill to 'terminate' all downstream channels
            sendPoisonPill()
            processor.terminate()
        }

    }

    final protected synchronized void sendPoisonPill() {

        taskConfig.outputs.each { name, channel ->
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
        taskConfig.outputs.keySet().each { name ->
            allFiles[ name ] = executor.collectResultFile(task, name)
        }

        return allFiles
    }

    synchronized protected bindOutputs( Map allOutputResources ) {

        log.trace "Binding results > task: ${name} - values: ${allOutputResources}"

        // -- bind each produced file to its own channel
        taskConfig.outputs.keySet().eachWithIndex { name, index ->

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

        // add the task specific entries
        if( taskConfig['env'] ) {
            taskConfig['env'].each { name, value -> result.put(name, value?.toString()) }
        }

        return result
    }



    /**
     * Execute the specified task shell script
     *
     * @param script The script string to be execute, e.g. a BASH script
     * @return {@code TaskDef}
     */
    final void launchTask( TaskRun task ) {
        executor.launchTask(task)
    }


    /**
     * Map used to delegate variable resolution to script scope
     */
    @Slf4j
    static class DelegateMap implements Map {

        @Delegate
        private Map<String,Object> local

        private BaseScript script

        DelegateMap(BaseScript script) {
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

