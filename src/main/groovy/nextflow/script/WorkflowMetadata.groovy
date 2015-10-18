/*
 * Copyright (c) 2013-2015, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2015, Paolo Di Tommaso and the respective authors.
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

package nextflow.script
import java.nio.file.Path
import java.nio.file.Paths

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.transform.ToString
import groovy.util.logging.Slf4j
import nextflow.Const
import nextflow.util.Duration
/**
 * Models workflow metadata properties and notification handler
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
@ToString(includeNames = true)
class WorkflowMetadata {

    /**
     * Project repository Git remote URL
     */
    String repository

    /**
     * Git commit ID of the executed repository
     */
    String commitId

    /**
     * Git branch/tag of the executed repository
     */
    String revision

    /**
     * Timestamp at workflow execution start
     */
    Date start

    /**
     * Timestamp at workflow execution complete
     */
    Date complete

    /**
     * Time elapsed to complete workflow execution
     */
    Duration duration

    /**
     * Docker image used to run workflow tasks. It report the docker image
     * defined in the pipeline configuration. When more than an image is used
     * it returns a {@link Map} object containing (process name, image name) pair entries.
     */
    def container

    /**
     * Command line as enter by the user to launch the workflow execution
     */
    String commandLine

    /**
     * Nextflow runtime info, it includes the following entries:
     * <li>version: runtime version number
     * <li>build: runtime build number
     * <li>timestamp: runtime compile timestamp
     */
    Map nextflow

    /**
     * Reports if the execution completed successfully
     */
    boolean success

    /**
     * Directory where workflow project is store on the computer
     */
    Path projectDir

    /**
     * Directory where workflow execution has been launched
     */
    Path launchDir

    /**
     * Workflow working directory
     */
    Path workDir

    /**
     * The exit status of the task that caused the workflow execution to fail
     */
    Integer exitStatus

    /**
     * Error message of the task that caused the workflow execution to fail
     */
    String errorMessage

    /**
     * Detailed error of the task that caused the workflow execution to fail
     */
    String errorReport

    final private ScriptRunner owner

    final private List<Closure> events = []

    /**
     * Initialise the available workflow properties
     *
     * @param owner An instance of {@link ScriptRunner}
     */
    WorkflowMetadata( ScriptRunner owner ) {
        this.owner = owner
        this.repository = owner.scriptFile.repository
        this.commitId = owner.scriptFile.commitId
        this.revision = owner.scriptFile.revision
        this.projectDir = owner.scriptFile.localPath
        this.start = new Date()
        this.container = owner.fetchContainers()
        this.commandLine = owner.commandLine
        this.nextflow = [version: Const.APP_VER, build: Const.APP_BUILDNUM, timestamp: Const.APP_TIMESTAMP_UTC]
        this.workDir = owner.session.workDir
        this.launchDir = Paths.get('.').complete()

        // check if there's a onComplete action in the config file
        registerConfigAction(owner.session.config.workflow as Map)
        owner.session.onShutdown { invokeOnComplete() }
    }

    /**
     * Implements the following idiom in the pipeline script:
     * <pre>
     *     window.onComplete {
     *         // do something
     *     }
     * </pre>
     *
     * @param action The action handler
     */
    void onComplete( Closure action ) {

        final clone = (Closure)action.clone()
        clone.delegate = owner.session.binding.variables
        clone.resolveStrategy = Closure.DELEGATE_ONLY

        events.add(clone)
    }

    /**
     * Implements the following idiom in the pipeline script:
     * <pre>
     *     window.onComplete = {
     *         // do something
     *     }
     * </pre>
     *
     * @param action The action handler
     */
    void setOnComplete( Closure action ) {
        events << action
    }

    /**
     * Register a onComplete handler defined in the nextflow config file
     *
     * @param workflowConfig
     */
    private void registerConfigAction( Map workflowConfig ) {
        if( !workflowConfig ) return
        if( workflowConfig.onComplete instanceof Closure ) {
            onComplete( (Closure)workflowConfig.onComplete )
        }
    }

    /**
     * Invoke the execution of the `onComplete` event handler(s)
     */
    @PackageScope
    void invokeOnComplete() {
        this.complete = new Date()
        this.duration = Duration.of( complete.time - start.time )
        this.success = !owner.session.aborted

        if( owner.session.fault ) {
            errorReport = owner.session.fault.report
            def task = owner.session.fault.task
            if( task ) {
                exitStatus = task.exitStatus != Integer.MAX_VALUE ? task.exitStatus : null
                def err = task.dumpStderr()
                if( !err ) err = task.dumpStdout()
                if( err ) errorMessage = err.join('\n')
            }
        }
        else {
            exitStatus = 0
        }

        events.each { Closure action ->
            try {
                action.call()
            }
            catch (Exception e) {
                log.error("Failed to invoke `workflow.onComplete` event handler", e)
            }
        }
    }

    /**
     * @return Render the workflow properties
     */
    String toString() {
        def result = new StringBuilder()
        result << 'repository: ' << repository
        result << ', projectDir: ' << projectDir
        result << ', commitId: ' << commitId
        result << ', revision: ' << revision
        result << ', startTime: ' << start
        result << ', endTime: ' << complete
        result << ', duration: ' << duration
        result << ', container: ' << container
        result << ', commandLine: ' << commandLine
        result << ', nextflow: ' << nextflow
        result << ', success: ' << success
        result << ', workDir: ' << workDir
        result << ', launchDir: ' << launchDir
        return result.toString()
    }

//    /**
//     * Define the {@code onComplete} event handle context. It allows to access metadata properties
//     * without having to specify the {@code workflow} prefix. If a property name does not exist in
//     * the {@link WorkflowMetadata} context it fallback on main script binding context
//     */
//    @CompileStatic
//    class EventHandlerContext {
//
//        final List<String> properties
//
//        EventHandlerContext() {
//            properties = WorkflowMetadata.this.metaClass.getProperties() *. name
//        }
//
//        Object getProperty(String name) {
//            try {
//                if( properties.contains(name) )
//                    return WorkflowMetadata.this.getProperty(name)
//                else
//                    throw new MissingPropertyException(name)
//            }
//            catch( MissingPropertyException e1 ) {
//                try {
//                    return owner.session.binding.getVariable(name)
//                }
//                catch ( MissingPropertyException e2 ) {
//                    return null
//                }
//            }
//        }
//
//    }
}