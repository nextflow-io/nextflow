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
import java.util.concurrent.ExecutorService
import java.util.concurrent.Executors
import java.util.concurrent.locks.Condition
import java.util.concurrent.locks.ReentrantLock

import com.hazelcast.config.Config
import com.hazelcast.core.Hazelcast
import com.hazelcast.core.HazelcastInstance
import com.hazelcast.core.IQueue
import com.hazelcast.core.ITopic
import groovy.util.logging.Slf4j
import nextflow.daemon.DaemonLauncher
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

    private IQueue<HzCmdCall> tasksQueue

    private ITopic<HzCmdResult> resultsTopic

    private ExecutorService executor

    private Map<UUID, HzRemoteSession> allSessions

    /**
     * Poor man singleton. The current daemon instance
     */
    static HzDaemon current

    HzDaemon() {  }

    @Override
    void launch(Map properties) {
        log.info "Launching Hazelcast daemon .. "
        System.setProperty('hazelcast.logging.type','slf4j')
        System.setProperty('hazelcast.system.log.enabled','true')
        def cfg = new Config()
        HzSerializerConfig.registerAll(cfg.getSerializationConfig())
        hazelcast = Hazelcast.newHazelcastInstance(cfg)

        tasksQueue = hazelcast.getQueue(TASK_SUBMITS_NAME)
        resultsTopic = hazelcast.getTopic(TASK_RESULTS_NAME)
        allSessions = hazelcast.getMap(SESSION_MAP)

        // -- set the current instance
        current = this

        // -- wait for tasks to be processed
        processTasks()
    }

    static ClassLoader getClassLoaderFor(UUID sessionId) {
        if( !current )
            throw new IllegalStateException('Missing daemon instance')

        def session = current.allSessions.get(sessionId)
        if( !session )
            throw new IllegalStateException("Missing remote session for ID: $sessionId")

        def loader = new GroovyClassLoader()
        session.classpath.each { URL url ->
            log.debug "Adding classpath url: $url"
            loader.addURL(url)
        }

        return loader
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

    protected void execute( HzCmdCall task ) {
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

    private publishResultBack( HzCmdCall cmd, result, Throwable error ) {
        log.trace "Publishing command result: $cmd"

        try {
            final obj = new HzCmdResult(cmd,result,error)
            resultsTopic.publish(obj)

        }
        catch( Throwable throwable ) {
            log.debug "Unable to publish command result: $cmd", throwable
            resultsTopic.publish( new HzCmdResult(cmd,result,throwable) )
        }

    }

}



