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

    Path localPath
    String repository
    String commitId
    String revision
    Date startTime
    Date endTime
    Duration duration
    def container
    String commandLine
    Map nextflow

    private ScriptRunner owner

    final private List<Closure> events = []

    WorkflowMetadata( ScriptRunner owner ) {
        this.owner = owner
        this.repository = owner.scriptFile.repository
        this.commitId = owner.scriptFile.commitId
        this.revision = owner.scriptFile.revision
        this.localPath = owner.scriptFile.localPath
        this.startTime = new Date()
        this.container = owner.fetchContainers()
        this.commandLine = owner.commandLine
        this.nextflow = [version: Const.APP_VER, build: Const.APP_BUILDNUM, timestamp: Const.APP_TIMESTAMP_UTC]

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
        events << action
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

    private void registerConfigAction( Map config ) {
        if( !config ) return
        if( config.onComplete instanceof Closure ) {
            events.add( (Closure)config.onComplete )
        }
    }

    /**
     * This method is required in order to allow the {@code onComplete} closure
     * defined in the {@code nextflow.config} file to access the {@code workflow}
     * object.
     *
     * @return The workflow object itself.
     */
    def getWorkflow() { return this }

    @PackageScope
    void invokeOnComplete() {
        this.endTime = new Date()
        this.duration = Duration.of( endTime.time - startTime.time )
        events.each { Closure action ->
            action.delegate = this
            action.resolveStrategy = Closure.DELEGATE_FIRST
            try {
                action.call()
            }
            catch (Exception e) {
                log.error("Cannot invoke `workflow.onComplete` action", e)
            }
        }
    }

    /**
     * @return Render the workflow properties
     */
    String toString() {
        def result = new StringBuilder()
        result << 'repository: ' << repository
        result << ', localPath: ' << localPath
        result << ', commitId: ' << commitId
        result << ', revision: ' << revision
        result << ', startTime: ' << startTime
        result << ', endTime: ' << endTime
        result << ', duration: ' << duration
        result << ', container: ' << container
        result << ', commandLine: ' << commandLine
        result << ', nextflow: ' << nextflow
        return result.toString()
    }
}