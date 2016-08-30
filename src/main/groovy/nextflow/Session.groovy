/*
 * Copyright (c) 2013-2016, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2016, Paolo Di Tommaso and the respective authors.
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

package nextflow
import static nextflow.Const.EXTRAE_TRACE_CLASS
import static nextflow.Const.S3_UPLOADER_CLASS

import java.lang.reflect.Method
import java.nio.file.Path
import java.nio.file.Paths
import java.util.concurrent.ConcurrentLinkedQueue
import java.util.concurrent.ExecutorService
import java.util.concurrent.Executors

import com.upplication.s3fs.S3OutputStream
import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import groovyx.gpars.GParsConfig
import groovyx.gpars.dataflow.operator.DataflowProcessor
import nextflow.dag.DAG
import nextflow.exception.AbortOperationException
import nextflow.exception.MissingLibraryException
import nextflow.file.FileHelper
import nextflow.processor.ErrorStrategy
import nextflow.processor.TaskDispatcher
import nextflow.processor.TaskFault
import nextflow.processor.TaskHandler
import nextflow.processor.TaskProcessor
import nextflow.script.ScriptBinding
import nextflow.trace.GraphObserver
import nextflow.trace.TimelineObserver
import nextflow.trace.TraceFileObserver
import nextflow.trace.TraceObserver
import nextflow.util.Barrier
import nextflow.util.ConfigHelper
import nextflow.util.Duration
import nextflow.util.NameGenerator

/**
 * Holds the information on the current execution
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class Session implements ISession {

    /**
     * Keep a list of all processor created
     */
    final List<DataflowProcessor> allProcessors = []

    /**
     * Dispatch tasks for executions
     */
    TaskDispatcher dispatcher

    /**
     * Holds the configuration object
     */
    Map config

    /**
     * Enable / disable tasks result caching
     */
    boolean cacheable

    /**
     * whenever it has been launched in resume mode
     */
    boolean resumeMode

    /**
     * The folder where tasks temporary files are stored
     */
    Path workDir

    /**
     * The folder where the main script is contained
     */
    Path baseDir

    /**
     * The pipeline script name (without parent path)
     */
    String scriptName

    /**
     * The class name used to compile the pipeline script
     */
    String scriptClassName

    /**
     * Mnemonic name of this run instance
     */
    String runName

    /**
     * Folder(s) containing libs and classes to be added to the classpath
     */
    List<Path> libDir

    private Path binDir

    /**
     * The unique identifier of this session
     */
    private UUID uniqueId

    private DAG dag

    private CacheDB cache

    private Barrier processesBarrier = new Barrier()

    private Barrier monitorsBarrier = new Barrier()

    private volatile boolean cancelled

    private volatile boolean aborted

    private volatile boolean terminated

    private volatile ExecutorService execService

    private volatile TaskFault fault

    private ScriptBinding binding

    private ClassLoader classLoader

    private Queue<Closure<Void>> shutdownCallbacks = new ConcurrentLinkedQueue<>()

    private int poolSize

    private Queue<TraceObserver> observers

    private Closure errorAction

    private boolean statsEnabled

    boolean getStatsEnabled() { statsEnabled }

    TaskFault getFault() { fault }

    /**
     * Creates a new session with an 'empty' (default) configuration
     */
    Session() {
        create(new ScriptBinding([:]))
    }

    /**
     * Create a new session given the {@link ScriptBinding} object
     *
     * @param binding
     */
    Session(ScriptBinding binding) {
        create(binding)
    }

    /**
     * Create a new session given the configuration specified
     *
     * @param config
     */
    Session(Map cfg) {
        final config = cfg instanceof ConfigObject ? cfg.toMap() : cfg
        create(new ScriptBinding(config))
    }

    /**
     * @return The current session {@link UUID}
     */
    UUID getUniqueId() { uniqueId }

    /**
     * @return The session max number of thread allowed
     */
    int getPoolSize() { poolSize }

    /**
     * @return The session {@link TaskDispatcher}
     */
    TaskDispatcher getDispatcher() { dispatcher }

    CacheDB getCache() { cache }

    /**
     * Creates a new session using the configuration properties provided
     *
     * @param binding
     */
    private void create( ScriptBinding binding ) {
        assert binding != null

        this.binding = binding
        this.config = binding.config

        // -- poor man session object dependency injection
        Global.setSession(this)
        Global.setConfig(config)

        // -- cacheable flag
        cacheable = config.cacheable

        // -- sets resumeMode and uniqueId
        if( config.resume ) {
            resumeMode = true
            uniqueId = UUID.fromString(config.resume as String)
        }
        else {
           uniqueId = UUID.randomUUID()
        }
        log.debug "Session uuid: $uniqueId"

        // -- set the run name
        this.runName = config.runName ?: NameGenerator.next()
        log.debug "Run name: $runName"

        // -- normalize taskConfig object
        if( config.process == null ) config.process = [:]
        if( config.env == null ) config.env = [:]

        if( !config.poolSize ) {
            final cpus = Runtime.getRuntime().availableProcessors()
            config.poolSize = cpus==1 ? 2 : cpus
        }

        // -- set the thread pool size
        this.poolSize = config.poolSize as int
        log.debug "Executor pool size: ${poolSize}"

        // -- create the task dispatcher instance
        this.dispatcher = new TaskDispatcher(this)

        // -- DGA object
        this.dag = new DAG(session:this)
    }

    /**
     * Initialize the session workDir, libDir, baseDir and scriptName variables
     */
    void init( Path scriptPath ) {

        this.workDir = ((config.workDir ?: 'work') as Path).complete()
        this.setLibDir( config.libDir as String )

        if(!workDir.mkdirs()) throw new AbortOperationException("Cannot create work-dir: $workDir -- Make sure you have write permissions or specify a different directory by using the `-w` command line option")
        log.debug "Work-dir: ${workDir} [${FileHelper.getPathFsType(workDir)}]"

        if( scriptPath ) {
            // the folder that contains the main script
            this.setBaseDir(scriptPath.parent)
            // set the script name attribute
            this.setScriptName(scriptPath.name)
        }

        this.observers = createObservers()
        this.statsEnabled = observers.size()>0

        cache = new CacheDB(uniqueId,runName).open()
    }

    /**
     * Given the `run` command line options creates the required {@link TraceObserver}s
     *
     * @param runOpts The {@code CmdRun} object holding the run command options
     * @return A list of {@link TraceObserver} objects or an empty list
     */
    @PackageScope
    Queue createObservers() {
        def result = new ConcurrentLinkedQueue()

        createTraceFileObserver(result)
        createTimelineObserver(result)
        createExtraeObserver(result)
        createDagObserver(result)

        return result
    }

    /**
     * create the Extrae trace observer
     */
    protected void createExtraeObserver(Collection<TraceObserver> result) {
        Boolean isEnabled = config.navigate('extrae.enabled') as Boolean
        if( isEnabled ) {
            try {
                result << (TraceObserver)Class.forName(EXTRAE_TRACE_CLASS).newInstance()
            }
            catch( Exception e ) {
                log.warn("Unable to load Extrae profiler",e)
            }
        }
    }

    /**
     * Create timeline report file observer
     */
    protected void createTimelineObserver(Collection<TraceObserver> result) {
        Boolean isEnabled = config.navigate('timeline.enabled') as Boolean
        if( isEnabled ) {
            String fileName = config.navigate('timeline.file')
            if( !fileName ) fileName = TimelineObserver.DEF_FILE_NAME
            def traceFile = (fileName as Path).complete()
            def observer = new TimelineObserver(traceFile)
            result << observer
        }
    }

    protected void createDagObserver(Collection<TraceObserver> result) {
        Boolean isEnabled = config.navigate('dag.enabled') as Boolean
        if( isEnabled ) {
            String fileName = config.navigate('dag.file')
            if( !fileName ) fileName = GraphObserver.DEF_FILE_NAME
            def traceFile = (fileName as Path).complete()
            def observer = new GraphObserver(traceFile)
            result << observer
        }
    }

    /*
     * create the execution trace observer
     */
    protected void createTraceFileObserver(Collection<TraceObserver> result) {
        Boolean isEnabled = config.navigate('trace.enabled') as Boolean
        if( isEnabled ) {
            String fileName = config.navigate('trace.file')
            if( !fileName ) fileName = TraceFileObserver.DEF_FILE_NAME
            def traceFile = (fileName as Path).complete()
            def observer = new TraceFileObserver(traceFile)
            config.navigate('trace.raw') { it -> observer.useRawNumbers(it == true) }
            config.navigate('trace.sep') { observer.separator = it }
            config.navigate('trace.fields') { observer.setFieldsAndFormats(it) }
            result << observer
    }
    }

    def Session start() {
        log.debug "Session start invoked"

        // register shut-down cleanup hooks
        Global.onShutdown { cleanUp() }
        // create tasks executor
        execService = Executors.newFixedThreadPool(poolSize)
        // signal start to tasks dispatcher
        dispatcher.start()
        // signal start to trace observers
        observers.each { trace -> trace.onFlowStart(this) }

        return this
    }

    ScriptBinding getBinding() { binding }

    ClassLoader getClassLoader() { classLoader }

    Session setClassLoader( ClassLoader loader ) {
        this.classLoader = loader
        return this
    }

    Barrier getBarrier() { monitorsBarrier }

    /**
     * The folder where script binaries file are located, by default the folder 'bin'
     * in the script base directory
     */
    Path getBinDir() {
        binDir
    }

    void setBaseDir( Path baseDir ) {
        this.baseDir = baseDir

        def path = baseDir.resolve('bin')
        if( path.exists() && path.isDirectory() ) {
            this.binDir = path
        }
        else {
            log.debug "Script base path does not exist or is not a directory: ${path}"
        }
    }

    def void setLibDir( String str ) {

        if( !str ) return

        def files = str.split( File.pathSeparator ).collect { String it -> Paths.get(it) }
        if( !files ) return

        libDir = []
        for( Path file : files ) {
            if( !file.exists() )
                throw new MissingLibraryException("Cannot find specified library: ${file.complete()}")

            libDir << file
        }
    }

    def List<Path> getLibDir() {
        if( libDir )
            return libDir

        libDir = []
        def localLib = baseDir ? baseDir.resolve('lib') : Paths.get('lib')
        if( localLib.exists() ) {
            log.debug "Using default localLib path: $localLib"
            libDir << localLib
        }
        return libDir
    }

    /**
     * Await the termination of all processors
     */
    void await() {
        log.debug "Session await"
        processesBarrier.awaitCompletion()
        log.debug "Session await > processes completed"
        terminated = true
        monitorsBarrier.awaitCompletion()
        log.debug "Session await > done"
    }

    void destroy() {
        log.trace "Session > destroying"

        if( !aborted ) {
            allProcessors *. join()
            log.trace "Session > after processors join"
        }

        cleanUp()
        log.trace "Session > after cleanup"

        execService.shutdown()
        execService = null
        log.trace "Session > executor shutdown"

        // -- close db
        cache?.close()

        // -- shutdown s3 uploader
        shutdownS3Uploader()

        log.debug "Session destroyed"
    }

    final protected void cleanUp() {

        log.trace "Shutdown: $shutdownCallbacks"
        while( shutdownCallbacks.size() ) {
            def hook = shutdownCallbacks.poll()
            try {
                if( hook )
                    hook.call()
            }
            catch( Exception e ) {
                log.debug "Failed to execute shutdown hook: $hook", e
            }
        }

        // -- invoke observers completion handlers
        while( observers.size() ) {
            def trace = observers.poll()
            try {
                if( trace )
                    trace.onFlowComplete()
            }
            catch( Exception e ) {
                log.debug "Failed to invoke observer completion handler: $trace", e
            }
        }

    }

    /**
     * Halt the pipeline execution choosing exiting immediately or completing current
     * pending task depending the chosen {@link ErrorStrategy}
     *
     * @param fault A {@link TaskFault} instance representing the error that caused the pipeline to stop
     */
    void fault(TaskFault fault, TaskHandler handler=null) {
        if( this.fault ) { return }
        this.fault = fault

        if( fault.strategy == ErrorStrategy.FINISH ) {
            cancel(handler)
        }
        else {
            abort(fault.error)
        }
    }

    /**
     * Cancel the pipeline execution waiting for the current running tasks to complete
     */
    @PackageScope
    void cancel(TaskHandler handler) {
        log.info "Execution cancelled -- Finishing pending tasks before exit"
        cancelled = true
        notifyError(handler)
        dispatcher.signal()
        processesBarrier.forceTermination()
        allProcessors *. terminate()
    }

    /**
     * Terminate the pipeline execution killing all running tasks
     *
     * @param cause A {@link Throwable} instance representing the execution that caused the pipeline execution to abort
     */
    void abort(Throwable cause = null) {
        if( aborted ) return
        log.debug "Session aborted -- Cause: ${cause}"
        aborted = true
        notifyError(null)
        dispatcher.signal()
        processesBarrier.forceTermination()
        monitorsBarrier.forceTermination()
        allProcessors *. terminate()
    }

    @PackageScope
    void forceTermination() {
        terminated = true
        processesBarrier.forceTermination()
        monitorsBarrier.forceTermination()
        allProcessors *. terminate()

        execService?.shutdownNow()
        GParsConfig.shutdown()
    }

    boolean isTerminated() { terminated }

    boolean isAborted() { aborted }

    boolean isCancelled() { cancelled }

    void processRegister(TaskProcessor process) {
        log.debug ">>> barrier register (process: ${process.name})"
        processesBarrier.register(process)
    }

    void processDeregister(TaskProcessor process) {
        log.debug "<<< barrier arrive (process: ${process.name})"
        processesBarrier.arrive(process)
    }

    DAG getDag() { this.dag }

    ExecutorService getExecService() { execService }

    /**
     * Register a shutdown hook to close services when the session terminates
     * @param Closure
     */
    void onShutdown( Closure shutdown ) {
        if( !shutdown )
            return

        shutdownCallbacks << shutdown
    }

    void notifyProcessCreate(TaskProcessor process) {
        for( TraceObserver it : observers ) {
            try {
                it.onProcessCreate(process)
            }
            catch( Exception e ) {
                log.debug(e.getMessage(), e)
            }
        }
    }


    /**
     * Notifies that a task has been submitted
     */
    void notifyTaskSubmit( TaskHandler handler ) {
        final task = handler.task
        log.info "[${task.hashLog}] ${task.runType.message} > ${task.name}"
        // -- save a record in the cache index
        cache.putIndexAsync(handler)

        for( TraceObserver it : observers ) {
            try {
                it.onProcessSubmit(handler)
            }
            catch( Exception e ) {
                log.debug(e.getMessage(), e)
            }
        }
    }

    /**
     * Notifies task start event
     */
    void notifyTaskStart( TaskHandler handler ) {
        for( TraceObserver it : observers ) {
            try {
                it.onProcessStart(handler)
            }
            catch( Exception e ) {
                log.debug(e.getMessage(), e)
            }
        }
    }

    /**
     * Notifies task termination event
     *
     * @param handler
     */
    void notifyTaskComplete( TaskHandler handler ) {
        // save the completed task in the cache DB
        cache.putTaskAsync(handler)

        // notify the event to the observers
        for( TraceObserver it : observers ) {
            try {
                it.onProcessComplete(handler)
            }
            catch( Exception e ) {
                log.debug(e.getMessage(), e)
            }
        }
    }


    void notifyTaskCached( TaskHandler handler ) {
        // -- save a record in the cache index
        cache.cacheTaskAsync(handler)

        for( TraceObserver it : observers ) {
            try {
                it.onProcessCached(handler)
            }
            catch( Exception e ) {
                log.error(e.getMessage(), e)
            }
        }
    }


    /**
     * Notify a task failure
     *
     * @param handler
     * @param e
     */
    void notifyError( TaskHandler handler ) {
        if( !errorAction )
            return

        try {
            errorAction.call( handler?.getTraceRecord() )
        }
        catch( Throwable e ) {
            log.debug(e.getMessage(), e)
        }
    }

    /**
     * Define the error event handler
     * @param action
     */
    void onError( Closure action ) {
        errorAction = action
    }


    @Memoized
    public getExecConfigProp( String execName, String name, Object defValue, Map env = null  ) {
        def result = ConfigHelper.getConfigProperty(config.executor, execName, name )
        if( result != null )
            return result

        // -- try to fallback sys env
        def key = "NXF_EXECUTOR_${name.toUpperCase().replaceAll(/\./,'_')}".toString()
        if( env == null ) env = System.getenv()
        return env.containsKey(key) ? env.get(key) : defValue
    }

    /**
     * Defines the number of tasks the executor will handle in a parallel manner
     *
     * @param execName The executor name
     * @param defValue The default value if setting is not defined in the configuration file
     * @return The value of tasks to handle in parallel
     */
    @Memoized
    public int getQueueSize( String execName, int defValue ) {
        getExecConfigProp(execName, 'queueSize', defValue) as int
    }

    /**
     * Determines how often a poll occurs to check for a process termination
     *
     * @param execName The executor name
     * @param defValue The default value if setting is not defined in the configuration file
     * @return A {@code Duration} object. Default '1 second'
     */
    @Memoized
    public Duration getPollInterval( String execName, Duration defValue = Duration.of('1sec') ) {
        getExecConfigProp( execName, 'pollInterval', defValue ) as Duration
    }

    /**
     *  Determines how long the executors waits before return an error status when a process is
     *  terminated but the exit file does not exist or it is empty. This setting is used only by grid executors
     *
     * @param execName The executor name
     * @param defValue The default value if setting is not defined in the configuration file
     * @return A {@code Duration} object. Default '90 second'
     */
    @Memoized
    public Duration getExitReadTimeout( String execName, Duration defValue = Duration.of('90sec') ) {
        getExecConfigProp( execName, 'exitReadTimeout', defValue ) as Duration
    }

    /**
     * Determines how often the executor status is written in the application log file
     *
     * @param execName The executor name
     * @param defValue The default value if setting is not defined in the configuration file
     * @return A {@code Duration} object. Default '5 minutes'
     */
    @Memoized
    public Duration getMonitorDumpInterval( String execName, Duration defValue = Duration.of('5min')) {
        getExecConfigProp(execName, 'dumpInterval', defValue) as Duration
    }

    /**
     * Determines how often the queue status is fetched from the cluster system. This setting is used only by grid executors
     *
     * @param execName The executor name
     * @param defValue  The default value if setting is not defined in the configuration file
     * @return A {@code Duration} object. Default '1 minute'
     */
    @Memoized
    public Duration getQueueStatInterval( String execName, Duration defValue = Duration.of('1min') ) {
        getExecConfigProp(execName, 'queueStatInterval', defValue) as Duration
    }

    static private void shutdownS3Uploader() {
        if( classWasLoaded(S3_UPLOADER_CLASS) ) {
            log.debug "AWS S3 uploader shutdown"
            S3OutputStream.shutdownExecutor()
        }
    }

    static private boolean classWasLoaded(String className) {
        Method find = ClassLoader.class.getDeclaredMethod("findLoadedClass", [String.class] as Class[] );
        find.setAccessible(true)
        return find.invoke(ClassLoader.getSystemClassLoader(), className)
    }


}
