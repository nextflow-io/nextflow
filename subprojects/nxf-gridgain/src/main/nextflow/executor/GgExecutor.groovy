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
import java.nio.file.Files
import java.nio.file.Path

import groovy.io.FileType
import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.exception.ProcessException
import nextflow.extension.FilesExtensions
import nextflow.processor.DelegateMap
import nextflow.processor.FileHolder
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
@CompileStatic
@ServiceName('gridgain')
@SupportedScriptTypes( [ScriptType.SCRIPTLET, ScriptType.GROOVY] )
class GgExecutor extends Executor {

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
            GgTaskHandler.createGroovyHandler(task, taskConfig,this)
        }
        else {
            GgTaskHandler.createScriptHandler(task, taskConfig, this)
        }

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


    /**
     * An adapter for GridGain compute task
     *
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
 * A task handler for GridGain  cluster
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class GgTaskHandler extends TaskHandler {

    private GgExecutor executor

    private ScriptType type

    private Path exitFile

    private Path outputFile


    /**
     * The result object for this task
     */
    private GridFuture future

    static GgTaskHandler createScriptHandler( TaskRun task, TaskConfig taskConfig, GgExecutor executor ) {
        def handler = new GgTaskHandler(task,taskConfig)
        handler.executor = executor
        handler.type = ScriptType.SCRIPTLET
        handler.exitFile = task.workDir.resolve(TaskRun.CMD_EXIT)
        handler.outputFile = task.workDir.resolve(TaskRun.CMD_OUTFILE)
        return handler
    }

    static GgTaskHandler createGroovyHandler( TaskRun task, TaskConfig taskConfig, GgExecutor executor ) {
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
            future = executor.execute( new GgBashTask(task) )
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

        if( isRunning() && (future.isCancelled() || (future.isDone() && (!exitFile || exitFile.lastModified()>0)))  ) {
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
 * Models a task executed remotely in a GridGain cluster node
 *
 * @param < T > The type of the value returned by the {@code #call} method
 */
@CompileStatic
abstract class GgBaseTask<T> implements GridCallable<T>, GridComputeJob {

    @GridLoggerResource
    private transient GridLogger log;

    /**
     * This field is used to transport the class attributes as a unique serialized byte array
     */
    private byte[] payload

    /**
     * Holds the class attributes in this map. Note: is defined as 'transient' because
     * the map content is serialized as a byte[] and saved to the {@code payload} field
     */
    private transient Map<String,Object> attrs = [:]

    /**
     * The local scratch dir where the task is actually executed in the remote node.
     * Note: is declared transient because it is valid only on the remote-side execution,
     * thus it do not need to be transported
     *
     */
    protected transient Path scratchDir

    /**
     * Create
     * @param task
     */
    protected GgBaseTask( TaskRun task ) {

        attrs.taskId = task.id
        attrs.name = task.name
        attrs.workDir = task.workDir
        attrs.targetDir = task.targetDir
        attrs.inputFiles = [:]
        attrs.outputFiles = []

        // -- The a mapping of input files and target names
        def allFiles = task.getInputFiles().values()
        for( List<FileHolder> entry : allFiles ) {
            if( entry ) for( FileHolder it : entry ) {
                attrs.inputFiles[ it.stagePath.name ] = it.storePath
            }
        }

        // -- the list of expected file names in the scratch dir
        attrs.outputFiles = task.getOutputFilesNames()

        payload = KryoHelper.serialize(attrs)
    }

    /** ONLY FOR TESTING PURPOSE */
    protected GgBaseTask() {}

    /**
     * @return The task unique ID
     */
    protected Object getTaskId() { attrs.taskId }

    /**
     * @return The task descriptive name (only for debugging)
     */
    protected String getName() { attrs.name }

    /**
     * @return The path where result files have to be copied
     */
    protected Path getTargetDir() { (Path)attrs.targetDir }

    /**
     * @return The task working directory i.e. the folder containing the scripts files, but
     * it is not the actual task execution directory
     */
    protected Path getWorkDir() { (Path)attrs.workDir }

    /**
     * @return The a mapping of input files and target names
     */
    Map<String,Path> getInputFiles() { (Map<String,Path>)attrs.inputFiles }

    /**
     * @return the list of expected file names in the scratch dir
     */
    List<String> getOutputFiles() { (List<String>)attrs.outputFiles }

    /**
     * Copies to the task input files to the execution folder, that is {@code scratchDir}
     * folder created when this method is invoked
     *
     */
    protected void stage() {

        if( attrs == null && payload )
            attrs = (Map<String,Object>)KryoHelper.deserialize(payload)

        // create a local scratch dir
        scratchDir = Files.createTempDirectory('nxf-task')

        if( !inputFiles )
            return

        // move the input files there
        for( Map.Entry<String,Path> entry : inputFiles.entrySet() ) {
            def target = scratchDir.resolve(entry.key)
            log?.debug "Task $name > staging path: '${entry.value}' to: '$target'"
            FilesExtensions.copyTo(entry.value, target)
        }
    }

    /**
     * Copy back the task output files from the execution directory in the local node storage
     * to the task {@code targetDir}
     */
    protected void unstage() {
        log?.debug "Unstaging file names: $outputFiles"

        if( !outputFiles )
            return

        // create a bash script that will copy the out file to the working directory
        if( !Files.exists(targetDir) )
            Files.createDirectories(targetDir)

        for( String name : outputFiles ) {
            try {
                copyToTargetDir(name)
            }
            catch( IOException e ) {
                log.error("Unable to copy result file: $name to target dir", e)
            }
        }

    }

    /**
     * Copy the file with the specified name from the task execution folder
     * to the {@code targetDir}
     *
     * @param fileName A file name relative to the {@code scratchDir}.
     *        It can contain the {@code *} and {@code ?} wildcards
     */
    protected void copyToTargetDir( String fileName ) {

        // TODO keep this aligned with "Executor#collectResultFile"

        String filePattern = fileName.replace("?", ".?").replace("*", ".*")

        // when there's not change in the pattern, try to find a single file
        if( filePattern == fileName ) {
            FilesExtensions.copyTo( scratchDir.resolve(fileName), targetDir )
        }
        else {
            scratchDir.eachFileMatch(FileType.ANY, ~/$filePattern/ ) { Path source ->
                FilesExtensions.copyTo( source, targetDir )
            }
        }
    }

    /**
     * Invoke the task execution. It calls the following methods in this sequence: {@code stage}, {@code execute0} and {@code unstage}
     *
     * @return The {@code execute0} result value
     * @throws ProcessException
     */
    @Override
    final T call() throws Exception {
        try {
            /*
             * stage the input files in the working are`
             */
            stage()

            /*
             * execute the task
             */
            final T result = execute0()

            /*
             * copy back the result files to the shared area
             */
            unstage()

            // return the exit status eventually
            return result
        }
        catch( Exception e ) {
            log.error("Cannot execute closure for task > $name", e)
            throw new ProcessException(e)
        }

    }

    /**
     * Just a synonym for {@code #call}
     *
     * @return The value returned by the task execution
     */
    final Object execute() {
        call()
    }

    /**
     * The actual task executor code provided by the extending subclass
     *
     * @return The value returned by the task execution
     */
    protected abstract T execute0()
}

/**
 * Execute a remote shell task into a remote GridGain cluster node
 */
@CompileStatic
class GgBashTask extends GgBaseTask<Integer>  {

    private static final long serialVersionUID = - 5552939711667527410L

    @GridLoggerResource
    private transient GridLogger log;

    private transient Process process

    /**
     * The command line to be executed
     */
    List shell

    String container

    Map environment

    Object stdin

    String script

    GgBashTask( TaskRun task ) {
        super(task)
        this.stdin = task.stdin
        this.container = task.container
        this.environment = task.processor.getProcessEnvironment()
        this.shell = task.processor.taskConfig.shell
        this.script = task.script
    }

    protected void unstage() {
        super.unstage()
        // copy the 'exit' file and 'output' file
        copyFromScratchToWorkDir(TaskRun.CMD_EXIT)
        copyFromScratchToWorkDir(TaskRun.CMD_OUTFILE)
    }


    private copyFromScratchToWorkDir( String name ) {
        try {
            Files.copy(scratchDir.resolve(name), workDir.resolve(name))
        }
        catch( Exception e ) {
            log.debug "Unable to copy file: '$name' from: '$scratchDir' to: '$workDir'"
        }
    }


    @Override
    protected Integer execute0() throws GridException {

        def wrapper = new BashWrapperBuilder(
                shell: shell,
                input: stdin,
                script: script,
                workDir: scratchDir,    // <-- important: the scratch folder is used as the 'workDir'
                targetDir: targetDir,
                container: container,
                environment: environment )
        shell.add( wrapper.build().toString() )

        log.debug "Running task > ${name}\n shell: ${shell.join(' ')}\n workdir: ${scratchDir.toFile()}"
        ProcessBuilder builder = new ProcessBuilder()
                .directory(scratchDir.toFile())
                .command(shell)
                .redirectErrorStream(true)

        // launch and wait
        process = builder.start()
        def result = process.waitFor()

        // make sure to destroy the process and close the streams
        try { process.destroy() }
        catch( Throwable e ) { }

        log.debug "Completed task > $name - exitStatus: $result"
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
        "${getClass().simpleName}[taskId: $taskId; name: $name; workDir: $workDir; scratchDir: ${scratchDir}]"
    }

}

/**
 * Execute a groovy closure task in a remote GridGain node
 */
@CompileStatic
class GgClosureTask extends GgBaseTask<GgResultData> {

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
     * The task closure serialized as a byte array
     */
    final byte[] codeObj

    /**
     * The task delegate context serialized as a byte array
     */
    final byte[] delegateObj


    GgClosureTask( TaskRun task, UUID sessionId ) {
        super(task)
        assert task
        this.sessionId = sessionId
        this.codeObj = SerializationUtils.serialize(task.code.dehydrate())
        this.delegateObj = (task.code.delegate as DelegateMap).dehydrate()
    }

    @Override
    protected GgResultData execute0() throws GridException {
        log.debug "Running closure for task > ${name}"

        def loader = provider.getClassLoaderFor(sessionId)
        def delegate = DelegateMap.rehydrate(delegateObj,loader)
        Closure closure = (Closure)InputStreamDeserializer.deserialize(codeObj,loader)
        Object result = closure.rehydrate(delegate, delegate.getScript(), delegate.getScript()).call()
        return new GgResultData(value: result, context: delegate?.getHolder())

    }


    @Override
    void cancel() {

    }

}

/**
 * Models the result of a remote closure task execution
 */
@CompileStatic
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

    Map getContext() {
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

