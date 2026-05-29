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

import nextflow.Session
import nextflow.SysEnv
import nextflow.util.Duration
import spock.lang.Specification
import spock.lang.Timeout
import test.TestHelper

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class LogsCheckpointTest extends Specification {

    def 'should configure default delay' () {
        given:
        def session = Mock(Session) {
            getWorkDir() >> TestHelper.createInMemTempDir()
            getConfig() >> [:]
        }
        and:
        def checkpoint = Spy(LogsCheckpoint)
        checkpoint.createHandler() >> Mock(LogsHandler)

        when:
        checkpoint.onFlowCreate(session)
        then:
        checkpoint.@interval == Duration.of('90s')

        cleanup:
        checkpoint.onFlowComplete()
    }

    def 'should configure delay via env var' () {
        given:
        SysEnv.push(TOWER_LOGS_CHECKPOINT_INTERVAL: '200s')
        def session = Mock(Session) {
            getWorkDir() >> TestHelper.createInMemTempDir()
            getConfig() >> [:]
        }
        and:
        def checkpoint = Spy(LogsCheckpoint)
        checkpoint.createHandler() >> Mock(LogsHandler)

        when:
        checkpoint.onFlowCreate(session)
        then:
        checkpoint.@interval == Duration.of('200s')

        cleanup:
        checkpoint.onFlowComplete()
        SysEnv.pop()
    }

    def 'should configure delay via config file' () {
        given:
        SysEnv.push(NXF_WORK: '/some/path', TOWER_LOGS_CHECKPOINT_INTERVAL: '200s')
        def session = Mock(Session) {
            getConfig()>>[tower:[logs:[checkpoint:[interval: '500s']]]]
            getWorkDir() >> TestHelper.createInMemTempDir()
        }
        and:
        def checkpoint = Spy(LogsCheckpoint)
        checkpoint.createHandler() >> Mock(LogsHandler)

        when:
        checkpoint.onFlowCreate(session)
        then:
        checkpoint.@interval == Duration.of('500s')

        cleanup:
        checkpoint.onFlowComplete()
        SysEnv.pop()
    }

    @Timeout(30)
    def 'checkpoint worker should be a daemon thread so it cannot keep the JVM alive' () {
        given:
        def session = Mock(Session) {
            getWorkDir() >> TestHelper.createInMemTempDir()
            getConfig() >> [:]
        }
        and:
        def checkpoint = Spy(LogsCheckpoint)
        checkpoint.createHandler() >> Mock(LogsHandler)

        when:
        checkpoint.onFlowCreate(session)
        then:
        checkpoint.@thread.isDaemon()

        cleanup:
        checkpoint.onFlowComplete()
    }

    @Timeout(30)
    def 'should checkpoint logs periodically' () {
        given:
        // count down on each save; await blocks until at least 3 saves happened
        def saves = new CountDownLatch(3)
        def handler = Mock(LogsHandler) {
            saveFiles() >> { saves.countDown() }
        }
        def session = Mock(Session) {
            getWorkDir() >> TestHelper.createInMemTempDir()
            getConfig() >> [tower:[logs:[checkpoint:[interval: '100ms']]]]
        }
        and:
        def checkpoint = Spy(LogsCheckpoint)
        checkpoint.createHandler() >> handler

        when:
        checkpoint.onFlowCreate(session)
        then:
        saves.await(10, TimeUnit.SECONDS)

        cleanup:
        checkpoint.onFlowComplete()
    }

    @Timeout(30)
    def 'should stop promptly when idle without waiting for the full interval' () {
        given:
        def handler = Mock(LogsHandler)
        def session = Mock(Session) {
            getWorkDir() >> TestHelper.createInMemTempDir()
            // a very long interval: stop must NOT wait for it
            getConfig() >> [tower:[logs:[checkpoint:[interval: '1h']]]]
        }
        and:
        def checkpoint = Spy(LogsCheckpoint)
        checkpoint.createHandler() >> handler
        and:
        checkpoint.onFlowCreate(session)
        // let the worker reach its wait state
        sleep(300)

        when:
        long t0 = System.currentTimeMillis()
        checkpoint.onFlowComplete()
        long elapsed = System.currentTimeMillis() - t0

        then:
        elapsed < 5_000
        !checkpoint.@thread.isAlive()
        // never reached a checkpoint within the 1h interval
        0 * handler.saveFiles()
    }

    @Timeout(30)
    def 'should not block shutdown when saveFiles is hung on a network call' () {
        given:
        def entered = new CountDownLatch(1)
        def release = new CountDownLatch(1) // never released -> saveFiles blocks forever
        def handler = Mock(LogsHandler) {
            saveFiles() >> { entered.countDown(); release.await() }
        }
        def session = Mock(Session) {
            getWorkDir() >> TestHelper.createInMemTempDir()
            getConfig() >> [tower:[logs:[checkpoint:[interval: '50ms', terminateTimeout: '500ms']]]]
        }
        and:
        def checkpoint = Spy(LogsCheckpoint)
        checkpoint.createHandler() >> handler
        and:
        checkpoint.onFlowCreate(session)
        // wait until the worker is actually stuck inside saveFiles
        assert entered.await(10, TimeUnit.SECONDS)

        when:
        long t0 = System.currentTimeMillis()
        checkpoint.onFlowComplete()
        long elapsed = System.currentTimeMillis() - t0

        then:
        // returned within ~terminateTimeout, NOT waiting for the never-ending saveFiles
        elapsed < 5_000
        // the stuck worker was abandoned rather than joined; it is a daemon so it
        // cannot keep the JVM alive even though it is still technically running
        checkpoint.@thread.isAlive()
        checkpoint.@thread.isDaemon()

        cleanup:
        release.countDown()
    }
}