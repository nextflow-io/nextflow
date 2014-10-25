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

package nextflow
import java.nio.file.Path
import java.nio.file.Paths
import java.util.concurrent.ExecutorService
import java.util.concurrent.Executors

import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.operator.DataflowProcessor
import jsr166y.Phaser
import nextflow.cli.CliOptions
import nextflow.cli.CmdRun
import nextflow.exception.MissingLibraryException
import nextflow.file.FileHelper
import nextflow.processor.TaskDispatcher
import nextflow.processor.TaskProcessor
import nextflow.trace.TraceFileObserver
import nextflow.trace.TraceObserver
import nextflow.util.ConfigHelper
import nextflow.util.Duration
/**
 * Holds the information on the current execution
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class Session {

    static final String EXTRAE_TRACE_CLASS = 'nextflow.extrae.ExtraeTraceObserver'

    /**
     * Keep a list of all processor created
     */
    final List<DataflowProcessor> allProcessors = []

    /**
     * Dispatch tasks for executions
     */
    final TaskDispatcher dispatcher

    /**
     * Holds the configuration object
     */
    def Map config

    /**
     * Enable / disable tasks result caching
     */
    def boolean cacheable

    /**
     * whenever it has been launched in resume mode
     */
    def boolean resumeMode

    /**
     * The folder where tasks temporary files are stored
     */
    def Path workDir = Paths.get('work').toAbsolutePath()

    /**
     * The folder where the main script is contained
     */
    def Path baseDir

    /**
     * The pipeline script name (without parent path)
     */
    def String scriptName

    /**
     * Folder(s) containing libs and classes to be added to the classpath
     */
    def List<Path> libDir

    /**
     * The unique identifier of this session
     */
    def final UUID uniqueId

    final private Phaser phaser = new Phaser()

    private volatile boolean aborted

    private volatile boolean terminated

    private volatile boolean delegateAbortToTaskMonitor

    private volatile ExecutorService execService

    final private List<Closure<Void>> shutdownCallbacks = []

    final int poolSize

    private List<TraceObserver> observers = []

    private boolean statsEnabled

    boolean getStatsEnabled() { statsEnabled }

    void delegateAbortToTaskMonitor( boolean value ) {
        this.delegateAbortToTaskMonitor = value
    }

    /**
     * Creates a new session with an 'empty' (default) configuration
     */
    def Session() {
        this([:])
    }


    /**
     * Creates a new session using the configuration properties provided
     *
     * @param config
     */
    def Session( Map config ) {
        assert config != null
        this.config = config instanceof ConfigObject ? config.toMap() : config

        // poor man singleton
        Global.setSession(this)
        Global.setConfig(config)

        // normalize taskConfig object
        if( config.process == null ) config.process = [:]
        if( config.env == null ) config.env = [:]

        // set unique session from the taskConfig object, or create a new one
        uniqueId = (config.session as Map)?.uniqueId ? UUID.fromString( (config.session as Map).uniqueId as String) : UUID.randomUUID()
        log.debug "Session uuid: $uniqueId"

        if( !config.poolSize ) {
            def cpus = Runtime.getRuntime().availableProcessors()
            config.poolSize = cpus >= 3 ? cpus-1 : 2
        }

        this.poolSize = config.poolSize as int
        System.setProperty('gpars.poolsize', config.poolSize as String)
        log.debug "Executor pool size: ${poolSize}"

        // create the task dispatcher instance
        dispatcher = new TaskDispatcher(this)
    }

    /**
     *  Initialize the session object by using the command line options provided by
     *  {@link CliOptions} object
     *
     */
    def void init( CmdRun runOpts, Path scriptFile ) {

        this.cacheable = runOpts.cacheable
        this.resumeMode = runOpts.resume != null

        // note -- make sure to use 'FileHelper.asPath' since it guarantee to handle correctly non-standard file system e.g. 'dxfs'
        this.workDir = FileHelper.asPath(runOpts.workDir).toAbsolutePath()
        this.setLibDir( runOpts.libPath )

        if( scriptFile ) {
            // the folder that contains the main script
            this.baseDir = scriptFile.parent
            // set the script name attribute
            this.scriptName = scriptFile.name
        }


        this.observers = createObservers( runOpts )
        this.statsEnabled = observers.size()>0
    }

    @PackageScope
    List createObservers( CmdRun runOpts ) {
        /*
         * create the execution trace observer
         */
        def result = []
        Boolean isEnabled = config.navigate('trace.enabled') as Boolean
        if( isEnabled || runOpts.withTrace ) {
            String fileName = runOpts.withTrace
            if( !fileName ) fileName = config.navigate('trace.file')
            if( !fileName ) fileName = TraceFileObserver.DEF_FILE_NAME
            def traceFile = Paths.get(fileName).complete()
            def observer = new TraceFileObserver(traceFile)
            config.navigate('trace.raw') { it -> observer.useRawNumbers(it == true) }
            config.navigate('trace.delim') { observer.delim = it }
            config.navigate('trace.fields') { observer.setFieldsAndFormats(it) }
            result << observer
        }

        /*
         * create the Extrae trace object
         */
        if( runOpts.withExtrae ) {
            try {
                result << (TraceObserver)Class.forName(EXTRAE_TRACE_CLASS).newInstance()
            }
            catch( Exception e ) {
                log.warn("Unable to load Extrae profiler ${Const.log_detail_tip_message}",e)
            }
        }

        return result
    }

    def Session start() {
        log.debug "Session start > phaser register (session)"

        /*
         * - register all of them in the dispatcher class
         * - register the onComplete event
         */
        for( TraceObserver trace : observers ) {
            log.debug "Registering observer: ${trace.class.name}"
            dispatcher.register(trace)
            onShutdown { trace.onFlowComplete() }
        }

        Global.onShutdown { cleanUp() }
        execService = Executors.newFixedThreadPool( poolSize )
        phaser.register()
        dispatcher.start()

        // signal start to trace observers
        observers.each { trace -> trace.onFlowStart(this) }

        return this
    }

    @PackageScope
    def getPhaser() { phaser }

    /**
     * The folder where script binaries file are located, by default the folder 'bin'
     * in the script base directory
     */
    @Memoized
    def Path getBinDir() {
        if( !baseDir ) {
            log.debug "Script base directory is null";
            return null
        }

        def path = baseDir.resolve('bin')
        if( !path.exists() || !path.isDirectory() ) {
            log.debug "Script base path does not exist or is not a directory: ${path}"
            return null
        }

        return path
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
        allProcessors *. join()
        terminated = true
        log.debug "<<< phaser deregister (session)"
        phaser.arriveAndAwaitAdvance()
        log.debug "Session await > done"
    }

    void destroy() {
        log.trace "Session destroying"
        if( execService ) execService.shutdown()
        cleanUp()
        log.debug "Session destroyed"
    }

    final synchronized protected void cleanUp() {
        log.trace "Shutdown: $shutdownCallbacks"
        List<Closure<Void>> all = new ArrayList<>(shutdownCallbacks)
        for( def hook : all ) {
            try {
                hook.call()
            }
            catch( Exception e ) {
                log.debug "Failed executing shutdown hook: $hook", e
            }
        }

        // -- after the first time remove all of them to avoid it's called twice
        shutdownCallbacks.clear()
    }

    void abort() {
        log.debug "Session abort -- terminating all processors"
        aborted = true
        allProcessors *. terminate()
        if( !delegateAbortToTaskMonitor )
            System.exit(ExitCode.SESSION_ABORTED)
    }

    boolean isTerminated() { terminated }

    boolean isAborted() { aborted }

    def int taskRegister(TaskProcessor process) {
        log.debug ">>> phaser register (process)"
        for( TraceObserver it : observers ) { it.onProcessCreate(process) }
        phaser.register()
    }

    def int taskDeregister(TaskProcessor process) {
        log.debug "<<< phaser deregister (process)"
        for( TraceObserver it : observers ) { it.onProcessDestroy(process) }
        phaser.arriveAndDeregister()
    }

    def ExecutorService getExecService() {
        execService
    }

    /**
     * Register a shutdown hook to close services when the session terminates
     * @param Closure
     */
    def void onShutdown( Closure shutdown ) {
        if( !shutdown )
            return

        shutdownCallbacks << shutdown
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


//    /**
//     * Create a table report of all executed or running tasks
//     *
//     * @return A string table formatted displaying the tasks information
//     */
//    String tasksReport() {
//
//        TableBuilder table = new TableBuilder()
//                .head('name')
//                .head('id')
//                .head('status')
//                .head('path')
//                .head('exit')
//
//        tasks.entries().each { Map.Entry<Processor, TaskDef> entry ->
//            table << entry.key.name
//            table << entry.value.id
//            table << entry.value.status
//            table << entry.value.workDirectory
//            table << entry.value.exitCode
//            table << table.closeRow()
//        }
//
//        table.toString()
//
//    }

}
