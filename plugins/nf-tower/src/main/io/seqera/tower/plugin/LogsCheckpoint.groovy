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

import java.util.concurrent.CountDownLatch
import java.util.concurrent.TimeUnit

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
 * log, report and timeline files.
 *
 * <h2>Concurrency design</h2>
 *
 * A single daemon worker thread loops forever, sleeping for {@code interval}
 * between checkpoints. The sleep is implemented as {@code stopLatch.await(interval)}
 * rather than {@code Thread.sleep}: when shutdown is requested we {@code countDown()}
 * the latch, which wakes the worker <em>immediately</em> without ever setting the
 * thread's interrupt flag. Avoiding the interrupt flag matters because cloud SDKs
 * (e.g. the AWS S3 client) observe it and abort in-flight uploads — that side effect
 * was the source of repeated bugs in earlier versions of this class.
 *
 * <h2>Why shutdown can never hang the head job</h2>
 *
 * {@link #onFlowComplete()} / {@link #onFlowError} run on the main shutdown thread.
 * They signal the worker and then {@code join} it for at most {@code terminateTimeout}.
 * If the worker is stuck in a hung network upload inside {@code saveFiles()} the join
 * times out and we <em>abandon</em> the worker. Because it is a daemon thread it cannot
 * keep the JVM alive, so the head job exits cleanly. The final, authoritative log upload
 * is performed elsewhere (see {@code CacheCommand}), so abandoning this best-effort
 * periodic uploader loses nothing.
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
    private final CountDownLatch stopLatch = new CountDownLatch(1)

    @Override
    void onFlowCreate(Session session) {
        this.session = session
        this.config = session.config
        this.handler = createHandler()
        this.interval = config.navigate('tower.logs.checkpoint.interval', defaultInterval()) as Duration
        this.terminateTimeout = config.navigate('tower.logs.checkpoint.terminateTimeout', defaultTerminateTimeout()) as Duration
        thread = Threads.start('tower-logs-checkpoint', this.&run)
    }

    protected LogsHandler createHandler() {
        new LogsHandler(session, SysEnv.get())
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

    /**
     * Signal the worker thread to stop and wait a bounded amount of time for it to
     * terminate. Never blocks shutdown indefinitely: if the worker is stuck in a hung
     * upload it is abandoned once {@code terminateTimeout} elapses.
     */
    protected void stop() {
        if( thread==null )
            return
        // wake the worker from its interval wait without touching the interrupt flag
        stopLatch.countDown()
        try {
            thread.join(terminateTimeout.toMillis())
        }
        catch (InterruptedException e) {
            Thread.currentThread().interrupt()
        }
        if( thread.isAlive() )
            log.warn "Logs checkpoint thread did not terminate within ${terminateTimeout} - abandoning it to allow the run to shut down"
    }

    protected void run() {
        log.debug "Starting logs checkpoint thread - interval: ${interval}"
        try {
            // await() returns true when a stop has been requested, false on timeout
            while( !await(interval) ) {
                handler.saveFiles()
            }
        }
        finally {
            log.debug "Terminating logs checkpoint thread"
        }
    }

    /**
     * Wait up to {@code interval} for a stop request.
     *
     * @return {@code true} if a stop was requested (the loop must terminate),
     *         {@code false} if the interval elapsed (time to checkpoint the logs).
     */
    protected boolean await(Duration interval) {
        try {
            return stopLatch.await(interval.toMillis(), TimeUnit.MILLISECONDS)
        }
        catch (InterruptedException e) {
            log.debug "Interrupted logs checkpoint thread"
            Thread.currentThread().interrupt()
            return true
        }
    }
}
