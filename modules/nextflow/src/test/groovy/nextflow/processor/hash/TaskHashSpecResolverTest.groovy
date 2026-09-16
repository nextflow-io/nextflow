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

import nextflow.SysEnv
import nextflow.processor.TaskRun
import spock.lang.Specification

class TaskHashSpecResolverTest extends Specification {

    def cleanup() {
        SysEnv.pop()
    }

    def 'defaults to std/v4 when nothing selects a spec'() {
        given:
        SysEnv.push([:])

        expect:
        TaskHashSpecResolver.defaultSpec().is(StdSpecs.STD_V4)
    }

    def 'NXF_TASK_HASH_VER selects a spec by id'() {
        given:
        SysEnv.push([NXF_TASK_HASH_VER: 'std/v2'])

        expect:
        TaskHashSpecResolver.defaultSpec().is(StdSpecs.STD_V2)
    }

    def 'an unknown id fails loudly rather than silently falling back'() {
        given:
        SysEnv.push([NXF_TASK_HASH_VER: 'std/nope'])

        when:
        TaskHashSpecResolver.defaultSpec()
        then:
        thrown(IllegalArgumentException)
    }

    def 'a plugin factory takes precedence over the default'() {
        given:
        SysEnv.push([:])
        def custom = new TaskHashSpec('plugin/v1', StdSpecs.STD_V4.bindings, StdSpecs.STD_V4.encoding)
        def factory = Mock(TaskHashSpecFactory)
        def task = Mock(TaskRun)

        when:
        def result = TaskHashSpecResolver.resolve(task, [factory])
        then:
        1 * factory.create(task) >> custom
        result.is(custom)
    }

    def 'a factory that abstains falls through to the next, then to the default'() {
        given:
        SysEnv.push([:])
        def abstaining = Mock(TaskHashSpecFactory)
        def task = Mock(TaskRun)

        when:
        def result = TaskHashSpecResolver.resolve(task, [abstaining])
        then:
        1 * abstaining.create(task) >> null
        result.is(StdSpecs.STD_V4)
    }
}
