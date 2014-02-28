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

import groovy.transform.EqualsAndHashCode
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.processor.DelegateMap
import nextflow.script.ScriptType
import org.apache.commons.lang.SerializationUtils
/**
 *  Hazelcast model constants
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
@EqualsAndHashCode
class HzRemoteSession implements Serializable {

    final UUID id

    final List<URL> classpath

    final String scriptClassName

    HzRemoteSession() { }

    HzRemoteSession(Session session) {
        id = session.getUniqueId()
        classpath = session.getClasspath()
        scriptClassName = session.getScriptClassName()
    }

    HzRemoteSession( UUID uuid, List<URL> classpath, String className ) {
        this.id = uuid
        this.classpath = classpath
        this.scriptClassName = className
    }
}


/**
 * Launches BASH task on a remote Hazelcast node
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@EqualsAndHashCode
class HzCmdCall implements Callable, Serializable {

    private static final long serialVersionUID = - 5552939711667527410L

    final UUID sender

    final taskId

    /* Note: use a File instead instead of a Path, because the latter is not serializable*/
    final File workDir

    /** The command line to be executed */
    final List commandLine

    /**
     * The context map when the command is a Groovy task
     */
    final Map context

    /**
     * The type of the script: Groovy native or Shell script
     */
    final ScriptType type

    /**
     * Store the closure as bye array
     */
    final byte[] codeBytes

    /**
     * The code to be execute the command is a Groovy task
     */
    private transient Closure codeObj

    /** ONLY FOR de-serialization purposes -- required by kryo serializer */
    protected HzCmdCall() {}

    HzCmdCall( UUID sender, taskId, Path workDir, List cmdLine ) {
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


    HzCmdCall( UUID sender, taskId, Path workDir, Closure code ) {
        assert sender
        assert taskId
        assert workDir
        assert code

        this.sender = sender
        this.taskId = taskId
        this.workDir = workDir.toFile()
        this.codeObj = code.dehydrate()
        this.type = ScriptType.GROOVY
        this.codeBytes = SerializationUtils.serialize(this.codeObj)

        if( code.delegate instanceof DelegateMap ) {
            this.context = (code.delegate as DelegateMap).getHolder()
        }
        else if( code.delegate instanceof Map) {
            this.context = (Map)code.delegate
        }

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
        if( codeObj == null ) {
            if( !codeBytes )
                throw new IllegalStateException('Missing closure byte-code')

            def loader = HzDaemon.getClassLoaderFor(sender)
            codeObj = InputStreamDeserializer.deserialize(codeBytes,loader)
        }

        codeObj.delegate = context
        codeObj.call()

    }

    String toString() {
        "${getClass().simpleName}[sender: $sender; taskId: $taskId; workDir: ${workDir}; cmdline: ${commandLine}]"
    }

}


/**
 * The result on a remote execution
 */
class HzCmdResult implements Serializable {

    private static final long serialVersionUID = - 6956540465417153122L ;

    /**
     * The session ID of the sender
     */
    final UUID sender

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
     * The update context (only for groovy closure task)
     */
    final Map context


    /**
     * Create a result starting from the command that originates it
     *
     * @param cmd
     * @param result
     * @param error
     */
    HzCmdResult( HzCmdCall cmd, result, Throwable error ) {
        this(cmd.sender, cmd.taskId, result, error, cmd.context)
    }

    HzCmdResult( UUID sender, taskId, result, Throwable error, Map context ) {
        this.sender = sender
        this.taskId = taskId
        this.value = result
        this.error = error
        this.context = context
    }

    /** ONLY FOR de-serialization purposes -- required by kryo serializer */
    protected HzCmdResult() {}

    String toString() {
        "${getClass().simpleName}[sender: $sender; taskId: $taskId; result: ${value}; error: ${error}]"
    }

    boolean equals(o) {
        if (this.is(o)) return true
        if (getClass() != o.class) return false

        HzCmdResult that = (HzCmdResult) o

        if (context != that.context) return false
        if (sender != that.sender) return false
        if (taskId != that.taskId) return false
        if (value != that.value) return false
        if (error?.class != that.error?.class) return false

        return true
    }

    int hashCode() {
        int result
        result = sender.hashCode()
        result = 31 * result + taskId.hashCode()
        result = 31 * result + (value != null ? value.hashCode() : 0)
        result = 31 * result + (error != null ? error.hashCode() : 0)
        result = 31 * result + (context != null ? context.hashCode() : 0)
        return result
    }
}

