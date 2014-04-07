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

package nextflow.executor
import java.nio.file.Path

import groovy.transform.EqualsAndHashCode
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.exception.ProcessException
import nextflow.processor.DelegateMap
import nextflow.processor.TaskConfig
import nextflow.processor.TaskHandler
import nextflow.processor.TaskMonitor
import nextflow.processor.TaskPollingMonitor
import nextflow.processor.TaskRun
import nextflow.script.ScriptType
import nextflow.util.Duration
import nextflow.util.InputStreamDeserializer
import nextflow.util.KryoHelper
import org.apache.commons.lang.SerializationUtils
import org.gridgain.grid.GridException
import org.gridgain.grid.GridFuture
import org.gridgain.grid.GridNode
import org.gridgain.grid.compute.GridComputeJob
import org.gridgain.grid.compute.GridComputeJobResult
import org.gridgain.grid.compute.GridComputeLoadBalancer
import org.gridgain.grid.compute.GridComputeTaskAdapter
import org.gridgain.grid.lang.GridCallable
import org.gridgain.grid.lang.GridInClosure
import org.gridgain.grid.logger.GridLogger
import org.gridgain.grid.resources.GridLoadBalancerResource
import org.gridgain.grid.resources.GridLoggerResource
import org.gridgain.grid.resources.GridUserResource
import org.jetbrains.annotations.Nullable

/**
 * A Nextflow executor based on GridGain services
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@SupportedScriptTypes( [ScriptType.SCRIPTLET, ScriptType.GROOVY] )
class GgExecutor extends AbstractExecutor {

    @PackageScope
    GgConnector connector

    /**
     * Initialize the executor by getting a reference to the Hazelcast connector
     */
    def void init() {
        super.init()
        connector = GgConnector.create(taskMonitor)
    }

    /**
     * Creates the task monitor for this executor
     * @return An instance of {@code TaskMonitor}
     */
    @Override
    protected TaskMonitor createTaskMonitor() {
        TaskPollingMonitor.create(session, name, Duration.of('1s'))
    }


    /**
     *  Creates an handler for the specified task
     */
    @Override
    TaskHandler createTaskHandler(TaskRun task) {

        if( task.type == ScriptType.GROOVY ) {
            return GgTaskHandler.createGroovyHandler(task, taskConfig,this)
        }

        /*
         * otherwise as a bash script
         */
        final bash = new BashWrapperBuilder(task)

        // staging input files
        bash.stagingScript = {
            final files = task.getInputFiles()
            final staging = stagingFilesScript(files)
            return staging
        }

        // unstage script
        bash.unstagingScript = {
            return unstageOutputFilesScript(task)
        }

        // create the wrapper script
        bash.build()

        GgTaskHandler.createScriptHandler(task, taskConfig, this)
    }


    TaskPollingMonitor getTaskMonitor() {
        (TaskPollingMonitor)super.getTaskMonitor()
    }


    GridFuture call( GridCallable command ) {
        connector.getCluster().compute().call(command)
    }

    GridFuture execute( GridComputeJob task ) {
        connector.getCluster().compute().execute( new GgTaskWrapper(task), null)
    }


    /*
     * link http://atlassian.gridgain.com/wiki/display/GG60/Load+Balancing
     */
    static class GgTaskWrapper extends GridComputeTaskAdapter  {

        // Inject load balancer.
        @GridLoadBalancerResource
        transient GridComputeLoadBalancer balancer

        private GridComputeJob theJob

        GgTaskWrapper( GridComputeJob job ) {
            this.theJob = job
        }

        @Override
        Map<? extends GridComputeJob, GridNode> map(List nodes, @Nullable Object arg) throws GridException {

            Map<GridComputeJob, GridNode> jobUnit = [:]
            jobUnit.put(theJob, balancer.getBalancedNode(theJob, null))
            return jobUnit
        }

        @Override
        Object reduce(List list) throws GridException {
            return list.get(0)
        }
    }

}


/**
 * A task handler for Hazelcast cluster
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class GgTaskHandler extends TaskHandler {

    private GgExecutor executor

    private Path exitFile

    private Path wrapperFile

    private Path outputFile

    private ScriptType type


    /**
     * The result object for this task
     */
    private GridFuture future

    static GgTaskHandler createScriptHandler( TaskRun task, TaskConfig taskConfig, GgExecutor executor ) {
        def handler = new GgTaskHandler(task,taskConfig)
        handler.executor = executor
        handler.exitFile = task.getCmdExitFile()
        handler.outputFile = task.getCmdOutputFile()
        handler.wrapperFile = task.getCmdWrapperFile()
        handler.type = ScriptType.SCRIPTLET
        return handler
    }

    static createGroovyHandler( TaskRun task, TaskConfig taskConfig, GgExecutor executor ) {
        def handler = new GgTaskHandler(task,taskConfig)
        handler.executor = executor
        handler.type = ScriptType.GROOVY
        return handler
    }

    private GgTaskHandler(TaskRun task, TaskConfig config) {
        super(task,config)
    }

    @Override
    void submit() {

        // submit to an hazelcast node for execution
        final sessionId = task.processor.session.uniqueId
        if( type == ScriptType.SCRIPTLET ) {
            final List cmdLine = new ArrayList(taskConfig.shell ?: 'bash' as List ) << wrapperFile.getName()
            future = executor.execute( new GgBashTask( task, cmdLine ) )
        }
        else {
            future = executor.execute( new GgClosureTask( task, sessionId ) )
        }

        future.listenAsync( { executor.getTaskMonitor().signalComplete(); } as GridInClosure )

        // mark as submitted -- transition to STARTED has to be managed by the scheduler
        status = TaskHandler.Status.SUBMITTED
        log.trace "Task $task > Submitted"
    }

    @Override
    boolean checkIfRunning() {
        if( isSubmitted() && future ) {
            log.trace "Task ${task} > RUNNING"
            status = TaskHandler.Status.RUNNING
            return true
        }

        return false
    }

    @Override
    boolean checkIfCompleted() {

        if( isRunning() && (future.isDone() || future.isCancelled()) && (!exitFile || exitFile.lastModified()>0) ) {
            status = TaskHandler.Status.COMPLETED

            final result = (GridComputeJobResult)future.get()
            if( result.getException() ) {
                task.error = result.getException()
                return true
            }

            if( result.isCancelled() ) {
                task.error = ProcessException.CANCELLED_ERROR
                return true
            }

            // -- the task output depend by the kind of the task executed
            if( isScriptlet() ) {
                task.stdout = outputFile
                task.exitStatus = result.getData() as Integer
            }
            else {
                def data = result.getData() as GgResultData
                task.stdout = data.value
                task.code.delegate = new DelegateMap( task.processor, data.context )
            }

            log.trace "Task ${task} > DONE"
            return true
        }

        return false
    }

    @Override
    void kill() {
        future?.cancel()
    }



    /**
     * @return Whenever is a shell script task
     */
    boolean isScriptlet() { type == ScriptType.SCRIPTLET }

    /**
     * @return Whenever is a groovy closure task
     */
    boolean isGroovy() { type == ScriptType.GROOVY }

}

/**
 * Execute a remote shell task
 */
class GgBashTask implements GridCallable<Integer>, GridComputeJob {

    private static final long serialVersionUID = - 5552939711667527410L

    @GridLoggerResource
    private transient GridLogger log;

    private transient Process process

    /**
     * Note: use a File instead instead of a Path, because the latter is not serializable
     */
    final File workDir

    /**
     * The command line to be executed
     */
    final List commandLine

    /**
     * The task name
     */
    final String name

    /**
     * The task id
     */
    final taskId


    GgBashTask( TaskRun task , List cmdLine ) {
        assert task
        assert cmdLine

        this.taskId = task.id
        this.name = task.name
        this.workDir = task.workDirectory.toFile()
        this.commandLine = cmdLine
    }


    @Override
    Integer call() throws Exception {
        execute() as Integer
    }


    @Override
    Object execute() throws GridException {
        log.debug "Launching script for task > ${name}"

        ProcessBuilder builder = new ProcessBuilder()
                .directory(workDir)
                .command(commandLine)
                .redirectErrorStream(true)

        process = builder.start()
        def result = process.waitFor()
        // make sure to destroy the process and close the streams
        try { process.destroy() }
        catch( Throwable e ) { }
        // return the exit value
        return result
    }

    @Override
    void cancel() {
        if( process ) {
            log.debug "Cancelling process for task > $name"
            process.destroy()
        }
        else {
            log.debug "No process to cancel for task > $name"
        }

    }


    String toString() {
        "${getClass().simpleName}[taskId: $taskId; name: $name; workDir: ${workDir}]"
    }

}

/**
 * Execute in a remote cluster note a groovy closure task
 */
class GgClosureTask implements GridComputeJob, GridCallable {

    private static final long serialVersionUID = 5515528753549263068L

    @GridLoggerResource
    private transient GridLogger log

    @GridUserResource
    private transient GgClassLoaderProvider provider

    /**
     * The client session identifier, it is required in order to access to
     * remote class-path
     */
    final UUID sessionId

    /**
     * The task unique ID
     */
    final taskId

    /**
     * The task descriptive name (only for debugging)
     */
    final String name

    /**
     * The task working directory
     */
    final File workDir

    /**
     * The task closure serialized as a byte array
     */
    final byte[] codeObj

    /**
     * The task delegate context serialized as a byte array
     */
    final byte[] delegateObj


    GgClosureTask( TaskRun task, UUID sessionId ) {
        assert task
        this.sessionId = sessionId
        this.taskId = task.id
        this.name = task.name
        this.workDir = task.workDirectory.toFile()
        this.codeObj = SerializationUtils.serialize(task.code.dehydrate())
        this.delegateObj = (task.code.delegate as DelegateMap).dehydrate()
    }

    @Override
    Object execute() throws GridException {
        log.debug "Running closure for task > ${name}"

        try {
            def loader = provider.getClassLoaderFor(sessionId)
            def delegate = DelegateMap.rehydrate(delegateObj,loader)
            Closure closure = InputStreamDeserializer.deserialize(codeObj,loader)
            Object result = closure.rehydrate(delegate, delegate.getScript(), delegate.getScript()).call()
            return new GgResultData(value: result, context: delegate?.getHolder())
        }
        catch( Exception e ) {
            log.error("Cannot execute closure for task > $name", e)
            throw new ProcessException(e)
        }
    }


    @Override
    void cancel() {

    }


    @Override
    GgResultData call() throws Exception {
        (GgResultData)execute()
    }

}

/**
 * Models the result of a remote closure task execution
 */
@EqualsAndHashCode
class GgResultData implements Serializable {

    private static final long serialVersionUID = - 7200781198107958188L ;

    /**
     * The closure returned value serialized as a byte array
     */
    private byte[] fValueObj

    /**
     * The closure execution context serialized as a byte array
     */
    private byte[] fContextObj

    transient Object value

    transient Map context

    Throwable error

    def getValue() {
        if( !value && fValueObj != null ) {
            value = KryoHelper.deserialize(fValueObj)
        }
        return value
    }

    void setValue( obj ) {
        this.value = obj
        this.fValueObj = KryoHelper.serialize(obj)
    }

    def getContext() {
        if( context == null && fContextObj != null ) {
            context = (Map)KryoHelper.deserialize(fContextObj)
        }
        return context
    }

    void setContext( Map ctx ) {
        this.context = ctx
        this.fContextObj = KryoHelper.serialize(ctx)
    }

}

