/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

package nextflow

import java.lang.reflect.Method
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths
import java.util.concurrent.ConcurrentLinkedQueue
import java.util.concurrent.ExecutorService
import java.util.concurrent.Executors

import com.google.common.hash.HashCode
import com.upplication.s3fs.S3OutputStream
import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import groovyx.gpars.GParsConfig
import groovyx.gpars.dataflow.operator.DataflowProcessor
import nextflow.config.Manifest
import nextflow.container.ContainerConfig
import nextflow.dag.DAG
import nextflow.exception.AbortOperationException
import nextflow.exception.AbortSignalException
import nextflow.exception.IllegalConfigException
import nextflow.exception.MissingLibraryException
import nextflow.file.FileHelper
import nextflow.processor.ErrorStrategy
import nextflow.processor.ProcessConfig
import nextflow.processor.TaskDispatcher
import nextflow.processor.TaskFault
import nextflow.processor.TaskHandler
import nextflow.processor.TaskProcessor
import nextflow.script.ScriptBinding
import nextflow.trace.GraphObserver
import nextflow.trace.ReportObserver
import nextflow.trace.StatsObserver
import nextflow.trace.AnsiLogObserver
import nextflow.trace.TimelineObserver
import nextflow.trace.TraceFileObserver
import nextflow.trace.TraceObserver
import nextflow.trace.TraceRecord
import nextflow.trace.WebLogObserver
import nextflow.trace.WorkflowStats
import nextflow.util.Barrier
import nextflow.util.ConfigHelper
import nextflow.util.Duration
import nextflow.util.HistoryFile
import nextflow.util.NameGenerator
import sun.misc.Signal
import sun.misc.SignalHandler
import static nextflow.Const.S3_UPLOADER_CLASS
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
    final Collection<DataflowProcessor> allOperators = new ConcurrentLinkedQueue<>()

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
     * Bucket work directory for cloud based executors
     */
    Path bucketDir

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

    /**
     * List files that concurrent on the session configuration
     */
    List<Path> configFiles

    /**
     * Local path where script generated classes are saved
     */
    private Path classesDir

    private Path binDir

    private Map<String,Path> binEntries = [:]

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

    private volatile Throwable error

    private ScriptBinding binding

    private ClassLoader classLoader

    private Queue<Closure<Void>> shutdownCallbacks = new ConcurrentLinkedQueue<>()

    private int poolSize

    private List<TraceObserver> observers

    private Closure errorAction

    private boolean statsEnabled

    private WorkflowStats workflowStats

    boolean getStatsEnabled() { statsEnabled }

    private boolean dumpHashes

    private List<String> dumpChannels

    boolean getDumpHashes() { dumpHashes }

    List<String> getDumpChannels() { dumpChannels }

    TaskFault getFault() { fault }

    Throwable getError() { error }

    WorkflowStats getWorkflowStats() { workflowStats }

    Path getClassesDir() { classesDir }

    boolean ansiLog

    private AnsiLogObserver ansiLogObserver

    AnsiLogObserver getAnsiLogObserver() { ansiLogObserver }

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
        this.dumpHashes = config.dumpHashes
        this.dumpChannels = (List<String>)config.dumpChannels

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
        log.debug "Work-dir: ${workDir.toUriString()} [${FileHelper.getPathFsType(workDir)}]"

        if( config.bucketDir ) {
            this.bucketDir = config.bucketDir as Path
            log.debug "Bucket-dir: ${bucketDir.toUriString()}"
        }

        if( scriptPath ) {
            // the folder that contains the main script
            this.setBaseDir(scriptPath.parent)
            // set the script name attribute
            this.setScriptName(scriptPath.name)
        }

        // set the byte-code target directory
        this.classesDir = FileHelper.createLocalDir()

        this.observers = createObservers()
        this.statsEnabled = observers.any { it.enableMetrics() }

        cache = new CacheDB(uniqueId,runName).open()
    }

    /**
     * Given the `run` command line options creates the required {@link TraceObserver}s
     *
     * @param runOpts The {@code CmdRun} object holding the run command options
     * @return A list of {@link TraceObserver} objects or an empty list
     */
    @PackageScope
    List<TraceObserver> createObservers() {
        def result = new ArrayList(10)

        createStatsObserver(result)     // keep as first, because following may depend on it
        createTraceFileObserver(result)
        createReportObserver(result)
        createTimelineObserver(result)
        createDagObserver(result)
        createWebLogObserver(result)
        createAnsiLogObserver(result)

        return result
    }

    protected void createAnsiLogObserver(Collection<TraceObserver> result) {
        if( ansiLog ) {
            this.ansiLogObserver = new AnsiLogObserver()
            result << ansiLogObserver
        }
    }

    /**
     * Create workflow message observer
     * @param result
     */
    protected void createWebLogObserver(Collection<TraceObserver> result) {
        Boolean isEnabled = config.navigate('weblog.enabled') as Boolean
        String url = config.navigate('weblog.url') as String
        if (isEnabled) {
            if ( !url ) url = WebLogObserver.DEF_URL
            def observer = new WebLogObserver(url)
            result << observer
        }
    }

    protected void createStatsObserver(Collection<TraceObserver> result) {
        final observer = new StatsObserver()
        this.workflowStats = observer.stats
        result << observer
    }


    /**
     * Create workflow report file observer
     */
    protected void createReportObserver(Collection<TraceObserver> result) {
        Boolean isEnabled = config.navigate('report.enabled') as Boolean
        if( isEnabled ) {
            String fileName = config.navigate('report.file')
            def maxTasks = config.navigate('report.maxTasks', ReportObserver.DEF_MAX_TASKS) as int
            if( !fileName ) fileName = ReportObserver.DEF_FILE_NAME
            def report = (fileName as Path).complete()
            def observer = new ReportObserver(report)
            observer.maxTasks = maxTasks

            result << observer
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

    /*
     * intercepts interruption signal i.e. CTRL+C
     * - on the first invoke session#abort
     * - on third force termination with System#exit
     */
    private void registerSignalHandlers() {

        int c = 0
        final ctrl_c = { Signal sig ->
            switch( ++c ) {
                case 1: abort(new AbortSignalException(sig)); println ''; break
                case 2: println "One more CTRL+C to force exit"; break
                default: log.info 'Adieu'; System.exit(1)
            }

        } as SignalHandler

        // -- abort session handler
        final abort_h = { Signal sig -> abort(new AbortSignalException(sig)) } as SignalHandler

        // -- register handlers
        Signal.handle( new Signal("INT"), ctrl_c)
        Signal.handle( new Signal("TERM"), abort_h)
        Signal.handle( new Signal("HUP"), abort_h)
    }

    /**
     * Dump the current dataflow network listing
     * the status of active processes and operators
     * for debugging purpose
     */
    String dumpNetworkStatus() {
        try {
            def msg = dag.dumpActiveNodes()
            msg ? "The following nodes are still active:\n" + msg : null
        }
        catch( Exception e ) {
            log.debug "Unexpected error dumping DGA status", e
            return null
        }
    }


    Session start() {
        log.debug "Session start invoked"

        // register shut-down cleanup hooks
        registerSignalHandlers()

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
    Path getBinDir() { binDir }

    Map<String,Path> getBinEntries() { binEntries ?: Collections.<String,Path>emptyMap() }

    void setBaseDir( Path baseDir ) {
        this.baseDir = baseDir

        def path = baseDir.resolve('bin')
        if( path.exists() && path.isDirectory() ) {
            this.binDir = path
            this.binEntries = findBinEntries(path)
        }
        else {
            log.debug "Script base path does not exist or is not a directory: ${path}"
        }
    }

    protected Map<String,Path> findBinEntries(Path path) {
        def result = [:]
        path
                .listFiles { file -> Files.isExecutable(file) }
                .each { Path file -> result.put(file.name,file)  }
        return result
    }

    void setLibDir( String str ) {

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

    List<Path> getLibDir() {
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

    @Memoized
    Manifest getManifest() {
        if( !config.manifest )
            return new Manifest()
        if( config.manifest instanceof Map )
            return new Manifest(config.manifest as Map)
        else {
            log.warn "Invalid config manifest definition [${this.getClass().getName()}]"
            return new Manifest()
        }
    }

    /**
     * Await the termination of all processors
     */
    void await() {
        log.debug "Session await"
        processesBarrier.awaitCompletion()
        log.debug "Session await > all process finished"
        terminated = true
        monitorsBarrier.awaitCompletion()
        log.debug "Session await > all barriers passed"
    }

    void destroy() {
        try {
            log.trace "Session > destroying"
            if( !aborted ) {
                allOperatorsJoin()
                log.trace "Session > after processors join"
            }

            shutdown0()
            log.trace "Session > after cleanup"

            execService.shutdown()
            execService = null
            log.trace "Session > executor shutdown"

            // -- close db
            cache?.close()

            // -- shutdown s3 uploader
            shutdownS3Uploader()

            // -- cleanup script classes dir
            classesDir.deleteDir()
        }
        finally {
            // -- update the history file
            if( HistoryFile.DEFAULT.exists() ) {
                HistoryFile.DEFAULT.update(runName,isSuccess())
            }
            log.trace "Session destroyed"
        }
    }

    final private allOperatorsJoin() {
        int attempts=0

        while( allOperators.size() ) {
            if( attempts++>0 )
                log.debug "This looks weird, attempt number $attempts to join pending operators"

            final itr = allOperators.iterator()
            while( itr.hasNext() ) {
                final op = itr.next()
                op.join()
                itr.remove()
            }
        }
    }


    final protected void shutdown0() {
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
        def copy = new ArrayList<TraceObserver>(observers)
        for( TraceObserver observer : copy  ) {
            try {
                if( observer )
                    observer.onFlowComplete()
            }
            catch( Exception e ) {
                log.debug "Failed to invoke observer completion handler: $observer", e
            }
        }

        // -- global
        Global.cleanUp()
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

        if( fault.task && fault.task.errorAction == ErrorStrategy.FINISH ) {
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
        allOperators *. terminate()
    }

    /**
     * Terminate the pipeline execution killing all running tasks
     *
     * @param cause A {@link Throwable} instance representing the execution that caused the pipeline execution to abort
     */
    void abort(Throwable cause = null) {
        if( aborted ) return
        log.debug "Session aborted -- Cause: ${cause?.message ?: cause ?: '-'}"
        aborted = true
        error = cause
        try {
            // log the dataflow network status
            def status = dumpNetworkStatus()
            if( status )
                log.debug(status)
            // force termination
            notifyError(null)
            dispatcher.signal()
            processesBarrier.forceTermination()
            monitorsBarrier.forceTermination()
            operatorsForceTermination()
        }
        catch( Throwable e ) {
            log.debug "Unexpected error while aborting execution", e
        }
    }

    private void operatorsForceTermination() {
        def operators = allOperators.toArray() as DataflowProcessor[]
        for( int i=0; i<operators.size(); i++ ) {
            operators[i].terminate()
        }
    }

    @PackageScope
    void forceTermination() {
        terminated = true
        processesBarrier.forceTermination()
        monitorsBarrier.forceTermination()
        allOperators *. terminate()

        execService?.shutdownNow()
        GParsConfig.shutdown()
    }

    boolean isTerminated() { terminated }

    boolean isAborted() { aborted }

    boolean isCancelled() { cancelled }

    boolean isSuccess() { !aborted && !cancelled }

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
     * Validate the config file
     *
     * @param processNames The list of process names defined in the pipeline script
     */
    void validateConfig(Collection<String> processNames) {
        def warns = validateConfig0(processNames)
        for( String str : warns )
            log.warn str
    }

    protected List<String> validateConfig0(Collection<String> processNames) {
        def result = []

        if( !(config.process instanceof Map) )
            return result

        // verifies that all process config names have a match with a defined process
        def keys = (config.process as Map).keySet()
        for(String key : keys) {
            String name = null
            if( key.startsWith('$') ) {
                name = key.substring(1)
            }
            else if( key.startsWith('withName:') ) {
                name = key.substring('withName:'.length())
            }
            if( name && !isValidProcessName(processNames, name, result) )
                break
        }

        return result
    }

    /**
     * Check that the specified name belongs to the list of existing process names
     *
     * @param selector The process name to check
     * @param processNames The list of processes declared in the workflow script
     * @param errorMessage A list of strings used to return the error message to the caller
     * @return {@code true} if the name specified belongs to the list of process names or {@code false} otherwise
     */
    protected boolean isValidProcessName(Collection<String> processNames, String selector, List<String> errorMessage)  {
        final matches = processNames.any { name -> ProcessConfig.matchesSelector(name, selector) }
        if( matches )
            return true

        def suggestion = processNames.closest(selector)
        def message = "There's no process matching config selector: $selector"
        if( suggestion )
            message += " -- Did you mean: ${suggestion.first()}?"
        errorMessage << message.toString()
        return false
    }
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
        for( int i=0; i<observers.size(); i++ ) {
            final observer = observers.get(i)
            try {
                observer.onProcessCreate(process)
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

        for( int i=0; i<observers.size(); i++ ) {
            final observer = observers.get(i)
            try {
                observer.onProcessSubmit(handler, handler.getTraceRecord())
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
        for( int i=0; i<observers.size(); i++ ) {
            final observer = observers.get(i)
            try {
                observer.onProcessStart(handler, handler.getTraceRecord())
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
        final trace = handler.getTraceRecord()
        cache.putTaskAsync(handler, trace)

        // notify the event to the observers
        for( int i=0; i<observers.size(); i++ ) {
            final observer = observers.get(i)
            try {
                observer.onProcessComplete(handler, trace)
            }
            catch( Exception e ) {
                log.debug(e.getMessage(), e)
            }
        }
    }


    void notifyTaskCached( TaskHandler handler ) {
        // -- save a record in the cache index
        cache.cacheTaskAsync(handler)

        for( int i=0; i<observers.size(); i++ ) {
            final observer = observers.get(i)
            try {
                observer.onProcessCached(handler, handler.getTraceRecord())
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

        for ( int i=0; i<observers?.size(); i++){
            try{
                final observer = observers.get(i)
                observer.onFlowError(handler, handler?.getTraceRecord())
            } catch ( Throwable e ) {
                log.debug(e.getMessage(), e)
            }
        }

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

    /**
     * Delete the workflow work directory from tasks temporary files
     */
    void cleanup() {
        if( !workDir || !config.cleanup )
            return

        if( aborted || cancelled || error )
            return

        CacheDB db = null
        try {
            log.trace "Cleaning-up workdir"
            db = new CacheDB(uniqueId, runName).openForRead()
            db.eachRecord { HashCode hash, TraceRecord record ->
                def deleted = db.removeTaskEntry(hash)
                if( deleted ) {
                    // delete folder
                    FileHelper.asPath(record.workDir).deleteDir()
                }
            }
            log.trace "Clean workdir complete"
        }
        catch( Exception e ) {
            log.warn("Failed to cleanup work dir: $workDir")
        }
        finally {
            db.close()
        }
    }

    /**
     * @return A {@link ContainerConfig} object representing the container engine configuration defined in config object
     */
    @Memoized
    ContainerConfig getContainerConfig() {

        def engines = new LinkedList<Map>()
        getContainerConfig0('docker', engines)
        getContainerConfig0('shifter', engines)
        getContainerConfig0('udocker', engines)
        getContainerConfig0('singularity', engines)

        def enabled = engines.findAll { it.enabled?.toString() == 'true' }
        if( enabled.size() > 1 ) {
            def names = enabled.collect { it.engine }
            throw new IllegalConfigException("Cannot enable more than one container engine -- Choose either one of: ${names.join(', ')}")
        }

        (enabled ? enabled.get(0) : ( engines ? engines.get(0) : [engine:'docker'] )) as ContainerConfig
    }


    private void getContainerConfig0(String engine, List<Map> drivers) {
        def config = this.config?.get(engine) as Map
        if( config ) {
            config.engine = engine
            drivers << config
        }
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

    @Memoized
    def getConfigAttribute(String name, defValue )  {
        def result = getMap0(getConfig(),name,name)
        if( result != null )
            return result

        def key = "NXF_${name.toUpperCase().replaceAll(/\./,'_')}".toString()
        def env = getSystemEnv()
        return (env.containsKey(key) ? env.get(key) : defValue)
    }

    private getMap0(Map map, String name, String fqn) {
        def p=name.indexOf('.')
        if( p == -1 )
            return map.get(name)
        else {
            def k=name.substring(0,p)
            def v=map.get(k)
            if( v == null )
                return null
            if( v instanceof Map )
                return getMap0(v,name.substring(p+1),fqn)
            throw new IllegalArgumentException("Not a valid config attribute: $fqn -- Missing element: $k")
        }
    }

    @Memoized
    protected Map<String,String> getSystemEnv() {
        new HashMap<String, String>(System.getenv())
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
    Duration getMonitorDumpInterval( String execName, Duration defValue = Duration.of('5min')) {
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
