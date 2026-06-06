/*
 * Copyright 2013-2026, Seqera Labs
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

package io.seqera.tower.plugin

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.SysEnv
import nextflow.trace.TraceObserverV2
import nextflow.trace.event.TaskEvent
import nextflow.util.Duration
import nextflow.util.Threads
/**
 * Implements a nextflow observer that periodically checkpoint
 * log, report and timeline files
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class LogsCheckpoint implements TraceObserverV2 {

    private Session session
    private Map config
    private Thread thread
    private Duration interval
    private Duration terminateTimeout
    private LogsHandler handler
    private volatile boolean stopped
    private final Object lock = new Object()

    @Override
    void onFlowCreate(Session session) {
        this.session = session
        this.config = session.config
        this.handler = new LogsHandler(session, SysEnv.get())
        this.interval = config.navigate('tower.logs.checkpoint.interval', defaultInterval()) as Duration
        this.terminateTimeout = config.navigate('tower.logs.checkpoint.terminateTimeout', defaultTerminateTimeout()) as Duration
        thread = Threads.start('tower-logs-checkpoint', this.&run)
    }

    private String defaultInterval() {
        SysEnv.get('TOWER_LOGS_CHECKPOINT_INTERVAL','90s')
    }

    private String defaultTerminateTimeout() {
        SysEnv.get('TOWER_LOGS_CHECKPOINT_TERMINATE_TIMEOUT','120s')
    }

    @Override
    void onFlowComplete() {
        stop()
    }

    @Override
    void onFlowError(TaskEvent event) {
        stop()
    }

    protected void stop() {
        if( thread == null )
            return
        synchronized(lock) {
            if( stopped )
                return
            stopped = true
            // wake up the checkpoint thread without relying on thread interruption
            lock.notifyAll()
        }
        // wait a bounded amount of time for the thread to terminate; if a saveFiles()
        // upload is genuinely hung the join times out and the daemon worker is abandoned
        // so it cannot keep the JVM alive and the run can shut down
        thread.join(terminateTimeout.toMillis())
        if( thread.isAlive() )
            log.warn "Logs checkpoint thread did not terminate within ${terminateTimeout} - abandoning it to allow the run to shut down"
    }

    protected void run() {
        log.debug "Starting logs checkpoint thread - interval: ${interval}"
        try {
            synchronized(lock) {
                while( !stopped ) {
                    // releases the lock and waits until the timeout elapses or stop() notifies
                    lock.wait(interval.toMillis())
                    if( stopped )
                        break
                    handler.saveFiles()
                }
            }
        }
        catch( InterruptedException e ) {
            log.debug "Interrupted logs checkpoint thread - cause: ${e.message}"
            Thread.currentThread().interrupt()
        }
        finally {
            log.debug "Terminating logs checkpoint thread"
        }
    }
}
