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

package nextflow.executor

import java.nio.file.Path
import java.util.concurrent.Callable

import groovy.transform.Canonical
import groovy.transform.EqualsAndHashCode
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.processor.DelegateMap
import nextflow.script.ScriptType

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
interface HzConst {

    final String TASK_RESULTS_NAME = 'results'

    final String TASK_SUBMITS_NAME = 'tasks'

    final String EXEC_SERVICE = 'default'

    final String SESSION_MAP = 'sessionMap'

}

/**
 * Keep track of the remote session classpath
 */
@Canonical
class HzRemoteSession implements Serializable {

    final UUID id

    final List<URL> classpath = []

    HzRemoteSession(Session session) {
        id = session.uniqueId
        classpath = session.getClasspath()
    }
}

/**
 * Launches BASH task on a remote Hazelcast node
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@EqualsAndHashCode
class HzBashCmd implements Callable, Serializable {

    private static final long serialVersionUID = - 5552939711667527410L

    final UUID sender

    final taskId

    /* Note: use a File instead instead of a Path, because the latter is not serializable*/
    final File workDir

    /** The command line to be executed */
    final List commandLine

    /**
     * The code to be execute the command is a Groovy task
     */
    final Closure code

    /**
     * The context map when the command is a Groovy task
     */
    final DelegateMap context

    /**
     * The type of the script: Groovy native or Shell script
     */
    final ScriptType type

    HzBashCmd( UUID sender, taskId, Path workDir, List cmdLine ) {
        assert sender
        assert taskId
        assert workDir
        assert cmdLine

        this.sender = sender
        this.taskId = taskId
        this.workDir = workDir.toFile()
        this.commandLine = cmdLine
        this.type = ScriptType.SCRIPTLET

    }


    HzBashCmd( UUID sender, taskId, Path workDir, Closure code ) {
        assert sender
        assert taskId
        assert workDir
        assert code

        this.sender = sender
        this.taskId = taskId
        this.workDir = workDir.toFile()
        this.context = (DelegateMap)code.delegate
        this.code = code.dehydrate()
        this.type = ScriptType.GROOVY

    }

    @Override
    Object call() throws Exception {
        if( type == ScriptType.SCRIPTLET )
            doBashCall()

        else
            doGroovyCall()
    }

    private Integer doBashCall( ) {
        log.trace "Launching script command > ${commandLine.join(' ')}"
        ProcessBuilder builder = new ProcessBuilder()
                .directory(workDir)
                .command(commandLine)
                .redirectErrorStream(true)

        def process = builder.start()
        return process.waitFor()
    }

    private doGroovyCall() {
        log.trace "Launching groovy command > ${taskId}"
        code.call()
    }

    String toString() {
        "${getClass().simpleName}[sender: $sender; taskId: $taskId; workDir: ${workDir}; cmdline: ${commandLine}]"
    }

}


/**
 * The result on a remote execution
 */
@EqualsAndHashCode
class HzBashResult implements Serializable {

    private static final long serialVersionUID = - 6956540465417153122L ;

    /**
     * The session ID of the sender
     */
    final UUID sender

    /**
     * The type of the command that originate this result
     */
    final ScriptType type

    /**
     * The unique task ID
     */
    final def taskId

    /**
     * The outcome of the executed command, either: the command exit status when the task is a shell
     * script of the closure return value when the task is a groovy command
     */
    final Object value

    /**
     * The task error if any
     */
    final Throwable error

    /**
     * @return Whenever is a shell script task
     */
    boolean isScriptlet() { type == ScriptType.SCRIPTLET }

    /**
     * @return Whenever is a groovy closure task
     */
    boolean isGroovy() { type == ScriptType.GROOVY }

    /**
     * Create a result starting from the command that originates it
     *
     * @param cmd
     * @param result
     * @param error
     */
    HzBashResult( HzBashCmd cmd, result, Throwable error ) {
        this(cmd.sender, cmd.type, cmd.taskId, result, error)
    }

    HzBashResult( UUID sender, ScriptType type, taskId, result, Throwable error ) {
        this.sender = sender
        this.type = type
        this.taskId = taskId
        this.value = result
        this.error = error
    }

    String toString() {
        "${getClass().simpleName}[sender: $sender; taskId: $taskId; type: $type; result: ${value}; error: ${error}]"
    }
}

