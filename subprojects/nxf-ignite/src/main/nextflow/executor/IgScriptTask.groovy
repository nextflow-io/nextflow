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

package nextflow.executor

import groovy.transform.CompileStatic
import nextflow.processor.TaskRun
import org.apache.ignite.IgniteException
import org.apache.ignite.IgniteLogger
import org.apache.ignite.resources.LoggerResource

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class IgScriptTask extends IgBaseTask<Integer>   {

    private static final long serialVersionUID = -1L

    private transient Process process

    @LoggerResource
    private transient IgniteLogger log

    /**
     * The command line to be executed
     */
    List<String> shell

    /**
     * The name of the container to run
     */
    String container

    /**
     * Whenever the process run an *executable* container
     */
    boolean executable

    Map environment

    Object stdin

    String script

    IgScriptTask( TaskRun task, UUID sessionId ) {
        super(task, sessionId)

        this.stdin = task.stdin
        this.container = task.container
        this.executable = task.isContainerExecutable()
        // note: create a copy of the process environment to avoid concurrent
        // process executions override each others
        this.environment = new HashMap( task.processor.getProcessEnvironment() )
        this.environment.putAll( task.getInputEnvironment() )
        this.shell = new ArrayList<>(task.config.getShell())
        this.script = task.script
    }

    @Override
    protected void stage() {
        // disable in-process staging
    }

    @Override
    protected void unstage() {
        // disable in-process un-staging
    }


    @Override
    protected Integer execute0() throws IgniteException {

        def session = getSessionFor(sessionId)
        if( log.isTraceEnabled() )
            log.trace "Session config: ${session.config}"


        def wrapper = new BashWrapperBuilder(
                shell: shell,
                input: stdin,
                script: script,
                scratch: scratch,
                workDir: workDir,
                targetDir: targetDir,
                container: container,
                environment: environment,
                copyStrategy: copyStrategy,
                dockerMount: localCacheDir,
                dockerConfig: session?.config?.docker,
                statsEnabled: session.statsEnabled,
                executable: executable
        )
        shell.add( wrapper.build().toString() )

        log.debug "Running task > name: ${name} - workdir: ${workDir}"
        ProcessBuilder builder = new ProcessBuilder()
                .directory(workDir.toFile())
                .command(shell)
                .redirectErrorStream(true)

        // launch and wait
        process = builder.start()
        def result = process.waitFor()

        // make sure to destroy the process and close the streams
        try { process.destroy() }
        catch( Throwable e ) { }

        log.debug "Completed task > $name - exitStatus: $result"
        // return the exit value
        return result
    }

    protected ScriptFileCopyStrategy getCopyStrategy() {
        new SimpleFileCopyStrategy(
                inputFiles: getInputFiles(),
                outputFiles: getOutputFiles(),
                targetDir: getTargetDir(),
                unstageStrategy: getUnstageStrategy()
        )
    }

    @Override
    void cancel() {
        if( process ) {
            log.debug "Cancelling process for task > $name"
            process.destroy()
        }
        else {
            log.debug "No process to cancel for task > $name"
        }
    }

    String toString() {
        "${getClass().simpleName}[taskId: $taskId; name: $name; workDir: $workDir]"
    }
}
