/*
 * Copyright 2013-2023, Seqera Labs
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
 *
 */

package io.seqera.tower.plugin

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.SysEnv
import nextflow.trace.TraceObserver
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
class LogsCheckpoint implements TraceObserver {

    private Session session
    private Map config
    private Thread thread
    private Duration interval
    private LogsHandler handler

    @Override
    void onFlowCreate(Session session) {
        this.session = session
        this.config = session.config
        this.handler = new LogsHandler(session, SysEnv.get())
        this.interval = config.navigate('tower.logs.checkpoint.interval', defaultInterval()) as Duration
    }

    private String defaultInterval() {
        SysEnv.get('TOWER_LOGS_CHECKPOINT_INTERVAL','90s')
    }

    @Override
    void onFlowBegin() {
        thread = Threads.start('tower-logs-checkpoint', this.&run)
    }

    @Override
    void onFlowComplete() {
        thread.interrupt()
        thread.join()
    }

    protected void run() {
        log.debug "Starting logs checkpoint thread - interval: ${interval}"
        try {
            while( !thread.isInterrupted() ) {
                // just wait the declared delay
                await(interval)
                // checkpoint the logs
                handler.saveFiles()
            }
        }
        finally {
            log.debug "Terminating logs checkpoint thread"
        }
    }

    protected void await(Duration interval) {
        try {
            Thread.sleep(interval.toMillis())
        }
        catch (InterruptedException e) {
            log.trace "Interrupted logs checkpoint thread"
            Thread.currentThread().interrupt()
        }
    }
}
