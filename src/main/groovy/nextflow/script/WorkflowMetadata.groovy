/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
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

import nextflow.config.Manifest

import java.nio.file.Path
import java.nio.file.Paths

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.transform.ToString
import groovy.util.logging.Slf4j
import nextflow.Const
import nextflow.config.ConfigBuilder
import nextflow.trace.WorkflowStats
import nextflow.util.Duration
import nextflow.util.VersionNumber
import org.codehaus.groovy.runtime.InvokerHelper
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
     * Run name
     */
    String runName

    /**
     * The script unique hash key
     */
    String scriptId

    /**
     * The main script file path
     */
    Path scriptFile

    /**
     * The main script name
     */
    String scriptName

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
     * User system home directory
     */
    Path homeDir

    /**
     * User system account name
     */
    String userName

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

    /**
     * The used configuration profile
     */
    String profile

    /**
     * The current sessionId
     */
    UUID sessionId

    /**
     * Returns ``true`` whenever the current instance is resumed from a previous execution
     */
    boolean resume

    /**
     * Which container engine was used to execute the workflow
     */
    String containerEngine

    /**
     * The list of files that concurred to create the config object
     */
    List<Path> configFiles

    /*
     * Workflow execution statistics
     */
    WorkflowStats stats

    /**
     * The workflow manifest
     */
    Manifest manifest

    final private ScriptRunner owner

    final private List<Closure> onCompleteActions = []

    final private List<Closure> onErrorActions = []

    /**
     * Initialise the available workflow properties
     *
     * @param owner An instance of {@link ScriptRunner}
     */
    WorkflowMetadata( ScriptRunner owner ) {
        this.owner = owner
        this.scriptId = owner.scriptFile.scriptId
        this.scriptFile = owner.scriptFile.main
        this.scriptName = owner.scriptFile.main?.fileName
        this.repository = owner.scriptFile.repository
        this.commitId = owner.scriptFile.commitId
        this.revision = owner.scriptFile.revision
        this.projectDir = owner.scriptFile.localPath
        this.start = new Date()
        this.container = owner.fetchContainers()
        this.commandLine = owner.commandLine
        this.nextflow = [version: new VersionNumber(Const.APP_VER), build: Const.APP_BUILDNUM, timestamp: Const.APP_TIMESTAMP_UTC]
        this.workDir = owner.session.workDir
        this.launchDir = Paths.get('.').complete()
        this.profile = owner.profile ?: ConfigBuilder.DEFAULT_PROFILE
        this.sessionId = owner.session.uniqueId
        this.resume = owner.session.resumeMode
        this.runName = owner.session.runName
        this.containerEngine = owner.session.containerConfig.with { isEnabled() ? getEngine() : null }
        this.configFiles = owner.session.configFiles?.collect { it.toAbsolutePath() }
        this.stats = owner.session.workflowStats
        this.userName = System.getProperty('user.name')
        this.homeDir = Paths.get(System.getProperty('user.home'))
        this.manifest = new Manifest((owner.session.config.manifest as Map) ?: Collections.emptyMap())

        // check if there's a onComplete action in the config file
        registerConfigAction(owner.session.config.workflow as Map)
        owner.session.onShutdown { invokeOnComplete() }
        owner.session.onError( this.&invokeOnError )
    }

    /**
     * Only for testing purpose -- do not use
     */
    @PackageScope
    WorkflowMetadata() {}

    /**
     * Implements the following idiom in the pipeline script:
     * <pre>
     *     workflow.onComplete {
     *         // do something
     *     }
     * </pre>
     *
     * @param action The action handler
     */
    void onComplete( Closure action ) {

        final clone = (Closure)action.clone()
        clone.delegate = owner.session.binding.variables
        clone.resolveStrategy = Closure.DELEGATE_FIRST

        onCompleteActions.add(clone)
    }

    /**
     * Implements the following idiom in the pipeline script:
     * <pre>
     *     workflow.onComplete = {
     *         // do something
     *     }
     * </pre>
     *
     * @param action The action handler
     */
    void setOnComplete( Closure action ) {
        onCompleteActions << action
    }

    /**
     * Implements the following idiom in the pipeline script:
     * <pre>
     *     workflow.onError {
     *         // do something
     *     }
     * </pre>
     * @param action
     */
    void onError( Closure action ) {

        final clone = (Closure)action.clone()
        clone.delegate = owner.session.binding.variables
        clone.resolveStrategy = Closure.DELEGATE_FIRST

        onErrorActions.add(clone)
    }

    /**
     * Dynamic getter for workflow metadata attributes
     *
     * @param field
     * @return The value associated to the specified field
     */
    def get(String field) {
        InvokerHelper.getProperty(this,field)
    }

    /**
     * Implements the following idiom in the pipeline script:
     * <pre>
     *     workflow.onError = {
     *         // do something
     *     }
     * </pre>
     *
     * @param action The action handler
     */
    void setOnError( Closure action ) {
        onErrorActions << action
    }

    /**
     * Register a onComplete handler defined in the nextflow config file
     *
     * @param workflowConfig
     */
    private void registerConfigAction( Map workflowConfig ) {
        if( !workflowConfig ) return
        // -- register `onComplete`
        if( workflowConfig.onComplete instanceof Closure ) {
            onComplete( (Closure)workflowConfig.onComplete )
        }
        // -- register `onError`
        if( workflowConfig.onError instanceof Closure ) {
            onError( (Closure)workflowConfig.onError )
        }
    }

    private void setErrorAttributes() {
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
        else if( owner.session.error ) {
            def msg = owner.session.error.message ?: owner.session.error.toString()
            errorMessage = msg
            errorReport = msg
        }
        else {
            exitStatus = 0
        }
    }

    /**
     * Invoke the execution of the `onComplete` event handler(s)
     */
    @PackageScope
    void invokeOnComplete() {
        this.complete = new Date()
        this.duration = Duration.of( complete.time - start.time )
        this.success = !(owner.session.aborted || owner.session.cancelled)

        setErrorAttributes()

        onCompleteActions.each { Closure action ->
            try {
                action.call()
            }
            catch (Exception e) {
                log.error("Failed to invoke `workflow.onComplete` event handler", e)
            }
        }

        // send email notification
        safeMailNotification()
    }

    void invokeOnError(trace) {
        this.success = false
        setErrorAttributes()
        onErrorActions.each { Closure action ->
            try {
                action.call(trace)
            }
            catch (Exception e) {
                log.error("Failed to invoke `workflow.onError` event handler", e)
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
        result << ', profile: ' << profile
        return result.toString()
    }

    /**
     * Tries to send the workflow completion email. Any exception is reported as a warning message.
     */
    protected void safeMailNotification() {
        try {
            def notifier = new WorkflowNotifier()
            notifier.workflow = this
            notifier.config = owner.session.config
            notifier.variables = owner.session.binding.variables
            notifier.sendNotification()
        }
        catch (Exception e) {
            log.warn "Failed to deliver notification email -- See the log file for details", e
        }
    }


}