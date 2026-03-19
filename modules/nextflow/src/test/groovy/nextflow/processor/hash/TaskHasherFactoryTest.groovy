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

package nextflow.processor.hash

import nextflow.Session
import nextflow.SysEnv
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TaskHasherFactoryTest extends Specification {

    def 'should create v1 hasher when configured'() {
        given:
        def session = Mock(Session) {
            getHashStrategy() >> TaskHasherFactory.Version.STD_V1
        }
        def processor = Mock(TaskProcessor) {
            getSession() >> session
        }
        def task = Mock(TaskRun) {
            getProcessor() >> processor
        }

        when:
        def hasher = TaskHasherFactory.create(task)

        then:
        hasher instanceof TaskHasherV1
        hasher.class == TaskHasherV1
    }

    def 'should create v2 hasher when configured'() {
        given:
        def session = Mock(Session) {
            getHashStrategy() >> TaskHasherFactory.Version.STD_V2
        }
        def processor = Mock(TaskProcessor) {
            getSession() >> session
        }
        def task = Mock(TaskRun) {
            getProcessor() >> processor
        }

        when:
        def hasher = TaskHasherFactory.create(task)

        then:
        hasher instanceof TaskHasherV2
        hasher.class == TaskHasherV2
    }

    def 'should resolve version from env var'() {
        given:
        SysEnv.push(['NXF_TASK_HASH_VER': 'std/v1'])

        expect:
        TaskHasherFactory.Version.DEFAULT() == TaskHasherFactory.Version.STD_V1

        cleanup:
        SysEnv.pop()
    }

    def 'should default to std/v2 when env var not set'() {
        given:
        SysEnv.push([:])

        expect:
        TaskHasherFactory.Version.DEFAULT() == TaskHasherFactory.Version.STD_V2

        cleanup:
        SysEnv.pop()
    }

    def 'should have correct string values'() {
        expect:
        TaskHasherFactory.Version.STD_V1.value == 'std/v1'
        TaskHasherFactory.Version.STD_V2.value == 'std/v2'
        TaskHasherFactory.Version.STD_V1.toString() == 'std/v1'
        TaskHasherFactory.Version.STD_V2.toString() == 'std/v2'
    }

    def 'should resolve version from string value'() {
        expect:
        TaskHasherFactory.Version.of('std/v1') == TaskHasherFactory.Version.STD_V1
        TaskHasherFactory.Version.of('std/v2') == TaskHasherFactory.Version.STD_V2
    }

    def 'should throw on unknown version string'() {
        when:
        TaskHasherFactory.Version.of('unknown')

        then:
        thrown(IllegalArgumentException)
    }
}
