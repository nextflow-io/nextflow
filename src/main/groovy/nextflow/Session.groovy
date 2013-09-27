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

import com.google.common.collect.LinkedHashMultimap
import com.google.common.collect.Multimap
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.Dataflow
import groovyx.gpars.dataflow.operator.DataflowProcessor
import groovyx.gpars.group.NonDaemonPGroup
import groovyx.gpars.group.PGroup
import groovyx.gpars.util.PoolUtils
import jsr166y.Phaser
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class Session {

    /**
     * Keep a list of all processor created
     */
    List<DataflowProcessor> allProcessors = []

    /**
     * Keep the list all executed tasks
     * note: LinkedHashMultimap preserves insertion order of entries, as well as the insertion order of keys, and the set of values associated with any one key.
     */
    Multimap<TaskProcessor, TaskRun> tasks = LinkedHashMultimap.create()

    /**
     * Holds the configuration object
     */
    def Map config

    /**
     * Enable / disable tasks result caching
     */
    def boolean cacheable

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
        if( !baseDir ) { log.debug "Script base directory is null"; return null }
        def path = new File(baseDir, 'bin')
        if( !path.exists() || !path.isDirectory() ) {
            return null
        }
        log.debug "Setting script bin dir: ${path}"
        return path
    }()

    /**
     * The unique identifier of this session
     */
    def final UUID uniqueId

    final private Phaser phaser = new Phaser()

    final private PGroup pgroup

    private boolean aborted


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

        // normalize taskConfig object
        if( config.task == null ) config.task = [:]
        if( config.env == null ) config.env = [:]

        // set unique session from the taskConfig object, or create a new one
        uniqueId = config.session?.uniqueId ? UUID.fromString( config.session.uniqueId.toString() ) : UUID.randomUUID()

        if( !config.poolSize ) {
            config.poolSize = PoolUtils.retrieveDefaultPoolSize()
        }

        log.debug "Executor pool size: ${config.poolSize}"

        // configure the dataflow thread group
        pgroup = new NonDaemonPGroup( config.poolSize as int )
        Dataflow.activeParallelGroup.set(pgroup)

        phaser.register()
    }

    /**
     * Await the termination of all processors
     */
    void await() {
        log.trace "Phaser await .. "
        phaser.arriveAndAwaitAdvance()
        log.debug "Phaser passed"
    }

    void terminate() {
        log.trace "Session join .. "
        allProcessors *. join()
        log.trace "Session joined"
        pgroup.shutdown()
        log.debug "Session terminated"
    }

    void abort() {
        log.debug "Session abort -- terminating all processors"
        aborted = true
        allProcessors *. terminate()
        System.exit( ExitCode.SESSION_ABORTED )
    }

    boolean isAborted() { aborted }

    def int taskRegister() {
        phaser.register()
    }

    def int taskDeregister() {
        phaser.arriveAndDeregister()
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
