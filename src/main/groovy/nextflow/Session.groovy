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

package nextflow
import java.nio.file.Path
import java.nio.file.Paths
import java.util.concurrent.ExecutorService
import java.util.concurrent.Executors

import groovy.transform.Memoized
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.Dataflow
import groovyx.gpars.dataflow.operator.DataflowProcessor
import groovyx.gpars.group.NonDaemonPGroup
import groovyx.gpars.group.PGroup
import groovyx.gpars.util.PoolUtils
import jsr166y.Phaser
import nextflow.processor.TaskDispatcher
import nextflow.util.Duration

/**
 * Holds the information on the current execution
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class Session {

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
     * The script name
     */
    def String scriptName = 'script1'

    /**
     * The folder where tasks temporary files are stored
     */
    def Path workDir = Paths.get('./work')

    /**
     * The folder where the main script is contained
     */
    def File baseDir

    /**
     * The folder where script binaries file are located, by default the folder 'bin'
     * in the script base directory
     */
    @Lazy
    File binDir = {
        if( !baseDir ) {
            log.debug "Script base directory is null";
            return null
        }
        def path = new File(baseDir, 'bin')
        if( !path.exists() || !path.isDirectory() ) {
            log.debug "Script base path does not exist or is not a directory: ${path}"
            return null
        }
        return path
    }()

    /**
     * The unique identifier of this session
     */
    def final UUID uniqueId

    final private Phaser phaser = new Phaser()

    final private PGroup pgroup

    private boolean aborted

    private boolean terminated

    private volatile ExecutorService execService

    /* Poor man singleton object */
    static Session currentInstance

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
        this.config = config

        // poor man singleton
        currentInstance = this

        // normalize taskConfig object
        if( config.process == null ) config.process = [:]
        if( config.env == null ) config.env = [:]

        // set unique session from the taskConfig object, or create a new one
        uniqueId = config.session?.uniqueId ? UUID.fromString( config.session.uniqueId.toString() ) : UUID.randomUUID()

        if( !config.poolSize ) {
            config.poolSize = PoolUtils.retrieveDefaultPoolSize()
        }

        // todo this settings seems to be unused -- check it
        if( !config.queueSize ) {
            config.queueSize = PoolUtils.retrieveDefaultPoolSize()
        }

        log.debug "Executor pool size: ${config.poolSize}"

        // configure the dataflow thread group
        pgroup = new NonDaemonPGroup( config.poolSize as int )
        Dataflow.activeParallelGroup.set(pgroup)

        dispatcher = new TaskDispatcher(this)
    }


    def Session start() {
        log.debug "Session start > phaser register (session) "
        phaser.register()
        dispatcher.start()
        return this
    }

    @PackageScope
    def getPhaser() { phaser }

    /**
     * Await the termination of all processors
     */
    void await() {
        log.debug "Session await > processors completion"
        allProcessors *. join()
        terminated = true
        log.debug "<<< phaser deregister (session)"
        phaser.arriveAndAwaitAdvance()
        log.debug "Session await > done"
    }

    void destroy() {
        log.trace "Session destroying"
        if( pgroup ) pgroup.shutdown()
        if( execService ) execService.shutdown()
        log.debug "Session destroyed"
    }

    void abort() {
        log.debug "Session abort -- terminating all processors"
        aborted = true
        allProcessors *. terminate()
        System.exit( ExitCode.SESSION_ABORTED )
    }

    boolean isTerminated() { terminated }

    boolean isAborted() { aborted }

    def int taskRegister() {
        log.debug ">>> phaser register (process)"
        phaser.register()
    }

    def int taskDeregister() {
        log.debug "<<< phaser deregister (process)"
        phaser.arriveAndDeregister()
    }

    def ExecutorService getExecService() {

        def local = execService
        if( local )
            return local

        synchronized(this) {
            if (execService == null) {
                execService = Executors.newCachedThreadPool()
            }
            return execService
        }
    }


    protected getExecConfigProp( String execName, String propName, Object defValue ) {

        def result = null

        // make sure that the *executor* is a map object
        // it could also be a plain string (when it specifies just the its name)
        if( config.executor instanceof Map ){
            if( execName && config.executor['$'+execName] instanceof Map ) {
                result = config.executor['$'+execName][propName]
            }

            if( result==null && config.executor[propName] ) {
                result = config.executor[propName]
            }
        }


        if( result==null ) {
            result = defValue
            log.trace "Undefined executor property: '$propName' -- fallback default value: $result"
        }

        return result

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
