/*
 * Copyright (c) 2013-2014, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2014, Paolo Di Tommaso and the respective authors.
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
import java.util.concurrent.Callable

import groovy.transform.Canonical
import groovy.transform.EqualsAndHashCode
import groovy.util.logging.Slf4j
import nextflow.Const
import nextflow.Session
import nextflow.exception.ProcessException
import nextflow.processor.DelegateMap
import nextflow.processor.TaskRun
import nextflow.script.ScriptType
import nextflow.util.InputStreamDeserializer
import org.apache.commons.lang.SerializationUtils
/**
 *  Hazelcast model constants
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
interface HzConst {

    final String TASK_SUBMITS_QUEUE = 'taskSubmits'

    final String ALL_NODES_MAP = 'allNodes'

    final String ALL_TASKS_MAP = 'allTasks'

    final String SESSIONS_MAP = 'sessionsMap'

    final String EXEC_SERVICE = 'executor'

    final String DEFAULT_GROUP_NAME = Const.APP_NAME

}

/**
 * Keep track of the remote session classpath
 */
@EqualsAndHashCode
class HzRemoteSession implements Serializable {

    private static final long serialVersionUID = - 3956143328835315200L ;

    final UUID id

    final List<URL> classpath

    /* this may be removed */
    final String scriptClassName

    HzRemoteSession() { }

    HzRemoteSession(Session session) {
        id = session.getUniqueId()
        classpath = session.getClasspath()
        scriptClassName = session.getScriptClassName()
    }

    /** only for test */
    protected HzRemoteSession( UUID uuid, List<URL> classpath  ) {
        this.id = uuid
        this.classpath = classpath
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

    final UUID sessionId

    final taskId

    /**
     * Note: use a File instead instead of a Path, because the latter is not serializable
     */
    final File workDir

    /**
     * The command line to be executed
     */
    final List commandLine


    /**
     * The type of the script: Groovy native or Shell script
     */
    final ScriptType type

    /**
     * Store the closure as byte array
     */
    final byte[] codeBytes

    /**
     * The delegate map serialized as an array of bytes
     */
    final byte[] delegateBytes

    /**
     * The code to be execute the command is a Groovy task
     */
    private transient Closure codeObj

    //private transient DelegateMap delegateMap


    /** ONLY FOR de-serialization purposes -- required by kryo serializer */
    protected HzCmdCall() {}

    HzCmdCall( UUID sessionId, TaskRun task , List cmdLine ) {
        assert task
        assert cmdLine

        this.sessionId = sessionId
        this.taskId = task.id
        this.workDir = task.workDirectory.toFile()
        this.commandLine = cmdLine
        this.type = ScriptType.SCRIPTLET

    }

    HzCmdCall( UUID sessionId, TaskRun task ) {
        assert task
        this.sessionId = sessionId
        this.taskId = task.id
        this.workDir = task.workDirectory.toFile()
        this.codeObj = task.code.dehydrate()
        //this.delegateMap = task.code.delegate as DelegateMap
        this.type = ScriptType.GROOVY
        // serialize to byte arrays
        //this.delegateBytes = this.delegateMap.dehydrate()
        this.codeBytes = SerializationUtils.serialize(this.codeObj)
    }

    //DelegateMap getDelegateMap() { delegateMap }

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

        // lookup for the class loader to be used to recreate the closure and the map holding the closure context e.g. variables
        def loader = HzDaemon.getClassLoaderFor(sessionId)

        if( delegateMap == null ) {
            delegateMap = DelegateMap.rehydrate(delegateBytes, loader)
            log.trace "Rehydrate delegate: $delegateMap"
        }

        if( codeObj == null ) {
            Closure closure = InputStreamDeserializer.deserialize(codeBytes,loader)
            // *rehydrate* the closure so that it can be executed
            codeObj = closure.rehydrate( delegateMap, delegateMap.getScript(), delegateMap.getScript() )
        }

        codeObj.call()

    }

    String toString() {
        "${getClass().simpleName}[session: $sessionId; taskId: $taskId; workDir: ${workDir}; cmdline: ${commandLine}]"
    }

}


/**
 * The result on a remote execution
 */
class HzCmdStatus implements Serializable {

    static final private long serialVersionUID = - 6956540465417153122L ;

    static enum Event { START, COMPLETE }

    /**
     * The type of this notification
     */
    final Event event

    /**
     * The ID of the session
     */
    final UUID sessionId

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
     * The member ID
     */
    final String memberId

    /**
     * just to make it unique
     */
    final long lastUpdate

    /**
     * Create a HzCmdStatus to notify a start event
     *
     * @param cmd
     * @param memberId
     * @return
     */
    static HzCmdStatus start( HzCmdCall cmd, String memberId ) {
        new HzCmdStatus( sessionId: cmd.sessionId, taskId: cmd.taskId, event: Event.START, memberId: memberId )
    }

    /**
     * Create a result starting from the command that originates it
     *
     * @param cmd
     * @param value
     * @param error
     */
    static HzCmdStatus result( HzCmdCall cmd, value, Throwable error ) {
        new HzCmdStatus( sessionId: cmd.sessionId, taskId: cmd.taskId, value: value, error: error, context: cmd.delegateMap?.getHolder() )
    }

    static HzCmdStatus error( UUID sessionId, taskId ) {
        new HzCmdStatus( sessionId: sessionId, taskId: taskId, error: ProcessException.UNKNOWN_ERROR )
    }


    /* NOTE: do not include the 'clock' attribute in the list of fields to copy */
    private static final FIELDS_TO_COPY = ['sessionId','taskId','value','error','context','event','memberId']

    final HzCmdStatus copyWith( Map params ) {
        final that = [:]
        for( String key : FIELDS_TO_COPY ) {
            that[key] = params.containsKey(key) ? params.get(key): this."$key"
        }

        new HzCmdStatus(that)
    }


    /**
     *
     * @param sessionId The client sessionId
     * @param taskId The {@code TaskRun#id}
     * @param value Task result when terminated
     * @param error Thrown exception (if any)
     * @param context Task evaluation context
     * @param event Kind of notification
     * @param memberId Cluster {@code Member} unique ID where the task is running
     */
    protected HzCmdStatus( Map params ) {
        this.sessionId = params.sessionId
        this.taskId = params.taskId
        this.value = params.value
        this.error = params.error
        this.context = params.context
        this.event = params.event ?: Event.COMPLETE
        this.memberId = params.memberId
        this.lastUpdate = params.clock ?: System.currentTimeMillis()
    }

    /** ONLY FOR de-serialization purposes -- required by kryo serializer */
    protected HzCmdStatus() {}

    String toString() {
        "${getClass().simpleName}[taskId: $taskId; event: $event; session: $sessionId; result: ${value}; error: ${error}; memberId: $memberId ]"
    }

    boolean equals(o) {
        if (this.is(o)) return true
        if (getClass() != o.class) return false

        HzCmdStatus that = (HzCmdStatus) o

        if (event != that.event) return false
        if (context != that.context) return false
        if (sessionId != that.sessionId) return false
        if (taskId != that.taskId) return false
        if (value != that.value) return false
        if (error?.class != that.error?.class) return false
        if (memberId != that.memberId ) return false
        if (lastUpdate != that.lastUpdate ) return false
        return true
    }

    int hashCode() {
        int result
        result = sessionId.hashCode()
        result = 31 * result + event.hashCode()
        result = 31 * result + taskId.hashCode()
        result = 31 * result + (this.value != null ? this.value.hashCode() : 0)
        result = 31 * result + (error != null ? error.hashCode() : 0)
        result = 31 * result + (context != null ? context.hashCode() : 0)
        result = 31 * result + (memberId != null ? memberId.hashCode() : 0)
        result = 31 * result + lastUpdate
        return result
    }

    boolean isStart() {
        event == Event.START
    }

    boolean isComplete() {
        event == Event.COMPLETE
    }
}


@Canonical
class HzNodeInfo implements Serializable {

    private static final long serialVersionUID = - 4449713933162037810L ;

    final int slots

}


@Canonical
class HzTaskKey implements Serializable {

    UUID sessionId;

    Object taskId

}


