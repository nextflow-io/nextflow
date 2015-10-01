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
import java.nio.file.Files

import groovy.transform.CompileStatic
import nextflow.processor.TaskRun
import org.apache.ignite.IgniteException
import org.apache.ignite.IgniteLogger
import org.apache.ignite.resources.LoggerResource
/**
 * Execute a remote shell task into a remote Ignite cluster node
 */
@CompileStatic
class IgBashTask extends IgBaseTask<Integer>  {

    private static final long serialVersionUID = - 5552939711667527410L

    @LoggerResource
    private transient IgniteLogger log

    private transient Process process

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

    IgBashTask( TaskRun task, UUID sessionId ) {
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

    protected void unstage() {
        super.unstage()
        // copy the 'exit' file and 'output' file
        copyFromScratchToWorkDir(TaskRun.CMD_EXIT)
        copyFromScratchToWorkDir(TaskRun.CMD_OUTFILE)
        copyFromScratchToWorkDir(TaskRun.CMD_ERRFILE, true)
        copyFromScratchToWorkDir(TaskRun.CMD_TRACE, true)
    }


    private void copyFromScratchToWorkDir( String name, boolean ignoreError = false ) {
        try {
            Files.copy(scratchDir.resolve(name), workDir.resolve(name))
        }
        catch( Exception e ) {
            if( !ignoreError )
                log.debug "Unable to copy file: '$name' from: '$scratchDir' to: '$workDir'"
        }
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
                workDir: scratchDir,    // <-- important: the scratch folder is used as the 'workDir'
                targetDir: targetDir,
                container: container,
                environment: environment,
                dockerMount: localCacheDir,
                dockerConfig: session?.config?.docker,
                statsEnabled: session.statsEnabled,
                executable: executable
        )
        shell.add( wrapper.build().toString() )

        log.debug "Running task > name: ${name} - workdir: ${scratchDir.toFile()}"
        ProcessBuilder builder = new ProcessBuilder()
                .directory(scratchDir.toFile())
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
        "${getClass().simpleName}[taskId: $taskId; name: $name; workDir: $workDir; scratchDir: ${scratchDir}]"
    }

}