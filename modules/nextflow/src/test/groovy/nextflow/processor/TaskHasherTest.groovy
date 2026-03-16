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

package nextflow.processor

import java.nio.file.Path

import nextflow.Session
import nextflow.SysEnv
import nextflow.script.ProcessConfig
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TaskHasherTest extends Specification {

    def 'should compute unique task hash' () {

        given:
        def session = Mock(Session) {
            getUniqueId() >> UUID.fromString('b69b6eeb-b332-4d2c-9957-c291b15f498c')
            getBinEntries() >> [:]
        }
        def processor = Mock(TaskProcessor) {
            getName() >> 'hello'
            getSession() >> session
            getConfig() >> Mock(ProcessConfig)
        }
        and:
        def task = Mock(TaskRun) {
            getSource() >> 'hello world'
            isContainerEnabled() >> false
            getContainer() >> null
            getConfig() >> Mock(TaskConfig)
            getProcessor() >> processor
        }
        and:
        def hasher = Spy(new TaskHasherV1(task))

        when:
        def uuid1 = hasher.compute()
        def uuid2 = hasher.compute()
        then:
        hasher.getTaskGlobalVars() >>> [
            [foo:'a', bar:'b'],
            [bar:'b', foo:'a']
        ]
        and: 'global vars should not affect task hash'
        uuid1 == uuid2
    }

    def 'should include referenced bin files in the task hash' () {

        given:
        def session = Mock(Session) {
            getBinEntries() >> [
                'foo.sh': Path.of('/some/path/foo.sh'),
                'bar.sh': Path.of('/some/path/bar.sh')
            ]
        }
        def processor = Mock(TaskProcessor) {
            getName() >> 'hello'
            getSession() >> session
        }
        def task = Mock(TaskRun) {
            getProcessor() >> processor
        }
        def hasher = new TaskHasherV1(task)

        when:
        def result = hasher.getTaskBinEntries('var=x foo.sh')
        then:
        result.size() == 1
        result.contains(Path.of('/some/path/foo.sh'))

        when:
        result = hasher.getTaskBinEntries('echo $(foo.sh); bar.sh')
        then:
        result.size() == 2
        result.contains(Path.of('/some/path/foo.sh'))
        result.contains(Path.of('/some/path/bar.sh'))
    }

    def 'should get task directive vars' () {
        given:
        def processor = Spy(TaskProcessor) {
            getConfig() >> Mock(ProcessConfig)
        }
        and:
        def config = new TaskConfig()
        config.cpus = 4
        config.ext.alpha = 'AAAA'
        config.ext.delta = { foo }
        config.ext.omega = "${-> bar}"
        config.context = [foo: 'DDDD', bar: 'OOOO']
        and:
        def task = Mock(TaskRun) {
            getConfig() >> config
            getProcessor() >> processor
            getVariableNames() >> { [ 'task.cpus', 'task.ext.alpha', 'task.ext.delta', 'task.ext.omega' ] as Set }
        }

        when:
        def result = new TaskHasherV1(task).getTaskExtensionDirectiveVars()
        then:
        result == [
            'task.ext.alpha': 'AAAA',
            'task.ext.delta': 'DDDD',
            'task.ext.omega': 'OOOO',
        ]
    }

    def 'should create v1 hasher when configured'() {
        given:
        def session = Mock(Session) {
            getHashStrategy() >> TaskHasherFactory.Version.V1
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

    def 'should create v2 hasher by default'() {
        given:
        def session = Mock(Session) {
            getHashStrategy() >> TaskHasherFactory.Version.V2
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
    }

    def 'should resolve version from env var'() {
        given:
        SysEnv.push(['NXF_TASK_HASH_VER': 'V1'])

        expect:
        TaskHasherFactory.Version.DEFAULT() == TaskHasherFactory.Version.V1

        cleanup:
        SysEnv.pop()
    }

    def 'should default to v2 when env var not set'() {
        given:
        SysEnv.push([:])

        expect:
        TaskHasherFactory.Version.DEFAULT() == TaskHasherFactory.Version.V2

        cleanup:
        SysEnv.pop()
    }

    def 'should compute hash entries for eval outputs'() {

        when:
        def result1 = TaskHasherV1.computeEvalOutputCommands([
            'nxf_out_eval_2': 'echo "value2"',
            'nxf_out_eval_1': 'echo "value1"',
            'nxf_out_eval_3': 'echo "value3"'
        ])

        def result2 = TaskHasherV1.computeEvalOutputCommands([
            'nxf_out_eval_3': 'echo "value3"',
            'nxf_out_eval_1': 'echo "value1"',
            'nxf_out_eval_2': 'echo "value2"'
        ])

        then:
        result1 == result2
        result1 == 'nxf_out_eval_1=echo "value1"\nnxf_out_eval_2=echo "value2"\nnxf_out_eval_3=echo "value3"'
    }
}
