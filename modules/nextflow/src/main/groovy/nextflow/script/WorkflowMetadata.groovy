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

package nextflow.script

import java.nio.file.Path
import java.nio.file.Paths
import java.time.OffsetDateTime

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.transform.ToString
import groovy.util.logging.Slf4j
import nextflow.NF
import nextflow.NextflowMeta
import nextflow.Session
import nextflow.config.ConfigBuilder
import nextflow.config.Manifest
import nextflow.trace.WorkflowStats
import nextflow.util.Duration
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
     *
     * Use OffDateTime instead of Date -- See https://stackoverflow.com/a/32443004/395921
     */
    OffsetDateTime start

    /**
     * Timestamp at workflow execution complete
     *
     * Use OffDateTime instead of Date -- See https://stackoverflow.com/a/32443004/395921
     */
    OffsetDateTime complete

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
    NextflowMeta nextflow

    /**
     * Reports if the execution completed successfully
     */
    boolean success

    /**
     * Directory where workflow project is store on the computer
     */
    Path projectDir

    /**
     * The name of the project when executed from a github repo or
     * the script name when executing a script file
     */
    String projectName

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

    private Session session

    final private List<Closure> onCompleteActions = []

    final private List<Closure> onErrorActions = []

    /**
     * Initialise the available workflow properties
     *
     * @param owner An instance of {@link ScriptRunner}
     */
    WorkflowMetadata( Session session, ScriptFile scriptFile ) {
        this.session = session
        this.scriptId = scriptFile?.scriptId
        this.scriptFile = scriptFile?.main
        this.scriptName = scriptFile?.main?.fileName
        this.repository = scriptFile?.repository
        this.commitId = scriptFile?.commitId
        this.revision = scriptFile?.revision
        this.projectDir = scriptFile?.localPath
        this.projectName = scriptFile?.projectName ?: scriptName
        this.start = OffsetDateTime.now()
        this.container = session.fetchContainers()
        this.commandLine = session.commandLine
        this.nextflow = NextflowMeta.instance
        this.workDir = session.workDir
        this.launchDir = Paths.get('.').complete()
        this.profile = session.profile ?:  ConfigBuilder.DEFAULT_PROFILE
        this.sessionId = session.uniqueId
        this.resume = session.resumeMode
        this.runName = session.runName
        this.containerEngine = session.containerConfig.with { isEnabled() ? getEngine() : null }
        this.configFiles = session.configFiles?.collect { it.toAbsolutePath() }
        this.stats = new WorkflowStats()
        this.userName = System.getProperty('user.name')
        this.homeDir = Paths.get(System.getProperty('user.home'))
        this.manifest = session.getManifest()

        // check if there's a onComplete action in the config file
        registerConfigAction(session.config.workflow as Map)
        session.onShutdown { invokeOnComplete() }
        session.onError( this.&invokeOnError )
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
        clone.delegate = NF.binding.variables
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
        clone.delegate = NF.binding.variables
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
        if( session.fault ) {
            errorReport = session.fault.report
            def task = session.fault.task
            if( task ) {
                exitStatus = task.exitStatus != Integer.MAX_VALUE ? task.exitStatus : null
                def err = task.dumpStderr()
                if( !err ) err = task.dumpStdout()
                if( err ) errorMessage = err.join('\n')
            }
        }
        else if( session.error ) {
            def msg = session.error.message ?: session.error.toString()
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
        this.complete = OffsetDateTime.now()
        this.duration = Duration.between( start, complete )
        this.success = !(session.aborted || session.cancelled)
        this.stats = getWorkflowStats()

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
        this.stats = getWorkflowStats()
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

    Map toMap() {
        final allProperties = this.metaClass.getProperties()
        final result = new LinkedHashMap(allProperties.size())
        for( MetaProperty property : allProperties ) {
            if( property.name == 'class' )
                continue
            try {
                result[property.name] = property.getProperty(this)
            }
            catch( GroovyRuntimeException e) {
                if( !e.message.startsWith('Cannot read write-only property') ) throw e
            }
        }
        return result
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
            notifier.config = session.config
            notifier.variables = NF.binding.variables
            notifier.sendNotification()
        }
        catch (Exception e) {
            log.warn "Failed to deliver notification email -- See the log file for details", e
        }
    }

    protected WorkflowStats getWorkflowStats() {
        session.statsObserver.getStats()
    }

}
