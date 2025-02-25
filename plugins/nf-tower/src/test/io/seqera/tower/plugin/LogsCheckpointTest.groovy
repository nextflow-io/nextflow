/*
 * Copyright 2013-2024, Seqera Labs
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

import nextflow.Session
import nextflow.SysEnv
import nextflow.util.Duration
import spock.lang.Specification
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
        def checkpoint = new LogsCheckpoint()

        when:
        checkpoint.onFlowCreate(session)
        then:
        checkpoint.@interval == Duration.of('90s')
    }

    def 'should configure delay via env var' () {
        given:
        SysEnv.push(TOWER_LOGS_CHECKPOINT_INTERVAL: '200s')
        def session = Mock(Session) {
            getWorkDir() >> TestHelper.createInMemTempDir()
            getConfig() >> [:]
        }
        and:
        def checkpoint = new LogsCheckpoint()

        when:
        checkpoint.onFlowCreate(session)
        then:
        checkpoint.@interval == Duration.of('200s')

        cleanup:
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
        def checkpoint = new LogsCheckpoint()

        when:
        checkpoint.onFlowCreate(session)
        then:
        checkpoint.@interval == Duration.of('500s')

        cleanup:
        SysEnv.pop()
    }
}
