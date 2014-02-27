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
import java.util.concurrent.ExecutorService
import java.util.concurrent.Executors
import java.util.concurrent.locks.Condition
import java.util.concurrent.locks.ReentrantLock

import com.hazelcast.config.Config
import com.hazelcast.core.Hazelcast
import com.hazelcast.core.HazelcastInstance
import com.hazelcast.core.IQueue
import com.hazelcast.core.ITopic
import groovy.transform.Canonical
import groovy.transform.EqualsAndHashCode
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.daemon.DaemonLauncher
import nextflow.processor.DelegateMap
import nextflow.script.ScriptType
/**
 * Run the Hazelcast daemon used to process user processes
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class HzDaemon implements HzConst, DaemonLauncher {

    private HazelcastInstance hazelcast

    private final int capacity = Runtime.getRuntime().availableProcessors()

    private ReentrantLock mutex

    private Condition notFull

    private volatile int count

    private IQueue<HzBashCmd> tasksQueue

    private ITopic<HzBashResult> resultsTopic

    private ExecutorService executor

    private Map<UUID, RemoteSession> allSessions

    HzDaemon() { }

    @Override
    void launch(Map properties) {
        log.info "Launching Hazelcast daemon .. "
        System.setProperty('hazelcast.logging.type','slf4j')
        System.setProperty('hazelcast.system.log.enabled','true')
        def cfg = new Config();
        hazelcast = Hazelcast.newHazelcastInstance(cfg)

        tasksQueue = hazelcast.getQueue(TASK_SUBMITS_NAME)
        resultsTopic = hazelcast.getTopic(TASK_RESULTS_NAME)
        allSessions = hazelcast.getMap(SESSION_MAP)

        // -- wait for tasks to be processed
        processTasks()
    }

    /**
     * Takes a task to be processed from the distributed queue an execute it
     */
    protected void processTasks() {

        // the shared queue
        mutex = new ReentrantLock()
        executor = Executors.newCachedThreadPool()

        Thread.start {
            count = 0

            // make sure to do not run more tasks than that the actual capacity
            while( true ) {
                mutex.withLock {
                    while( count >= capacity ) { notFull.await() }
                    count++
                }

                execute(tasksQueue.take())
            }

        }
    }

    protected void execute( HzBashCmd task ) {
        log.trace "Executing task command > $task"

        executor.submit( {

            def result = null
            Throwable error = null
            try {
                result = task.call()
            }
            catch( Throwable ex ) {
                error = ex
            }
            finally {
                publishResultBack(task,result,error)
                mutex.withLock {
                    count--;
                    notFull.signal()
                }
            }


        } as Runnable )

    }

    private publishResultBack( HzBashCmd cmd, result, Throwable error ) {
        log.trace "Publishing command result: $cmd"

        try {
            final obj = new HzBashResult(cmd,result,error)
            resultsTopic.publish(obj)

        }
        catch( Throwable throwable ) {
            log.debug "Unable to publish command result: $cmd", throwable
            resultsTopic.publish( new HzBashResult(cmd,result,throwable) )
        }

    }


}

/**
 * Keep track of the remote session classpath
 */
@Canonical
class RemoteSession implements Serializable {

    final UUID id

    final List<URL> classpath = []

    RemoteSession(Session session) {
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
     * @return Wheneven is a groovy closure task
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


