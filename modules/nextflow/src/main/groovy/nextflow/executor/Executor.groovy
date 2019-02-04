/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

package nextflow.executor

import java.nio.file.Path

import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.file.FileHelper
import nextflow.processor.TaskHandler
import nextflow.processor.TaskMonitor
import nextflow.processor.TaskRun
import nextflow.script.ScriptType
/**
 * Declares methods have to be implemented by a generic
 * execution strategy
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@SupportedScriptTypes( [ScriptType.SCRIPTLET] )
abstract class Executor {

    /**
     * The current session object
     */
    Session session

    /**
     * The executor simple name
     */
    String name

    /**
     * The queue holder that keep track of all tasks for this executor.
     */
    private TaskMonitor monitor

    private static final Map<Class,Boolean> registerFlag = [:]

    /**
     * Template method executed once for executor class
     */
    void register() { }

    protected String getDisplayName() { name }

    /**
     * Allows to post-initialize the executor
     */
    void init() {
        log.debug "Initializing executor: $name"

        // -- skip if already assigned, this is only for testing purpose
        if( monitor )
            return

        // -- get the reference to the monitor class for this executor
        monitor = session.dispatcher.getOrCreateMonitor(this.class) {
            log.info "[warm up] executor > ${getDisplayName()}"
            createTaskMonitor()
        }

        // call the register template method
        if( !registerFlag[this.class] ) {
            log.debug "Invoke register for executor: ${name}"
            register()
            registerFlag[this.class] = true
        }
    }

    /**
     * The path where scratch data is written for the current executor.
     *
     * @return The executor base work directory
     */
    Path getWorkDir() {
       session.getWorkDir()
    }

    /**
     * Temporary work directory relative to the executor work directory
     *
     * @return The temporary directory path
     */
    Path getTempDir( String name = null, boolean create = true ) {
        def path = FileHelper.createTempFolder(getWorkDir())
        if( name )
            path = path.resolve(name)

        if( !path.exists() && create && !path.mkdirs() )
            throw new IOException("Unable to create folder: $path -- Check file system permission" )

        return path
    }

    /**
     * @return Create a new instance of the {@code TaskQueueHolder} component
     */
    abstract protected TaskMonitor createTaskMonitor()

    /**
     * @return A reference to the current {@code #queueHolder} object
     */
    @PackageScope
    TaskMonitor getTaskMonitor()  { monitor }

    /**
     * @return Create a new {@code TaskHandler} to manage the scheduling
     * actions for this task
     */
    abstract TaskHandler createTaskHandler(TaskRun task)

    /**
     * @return {@code true} whenever the containerization is managed by the executor itself
     */
    boolean isContainerNative() {
        return false
    }

}
