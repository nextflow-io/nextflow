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
 */

package nextflow.processor

import java.nio.file.Paths

import com.google.common.hash.HashCode
import nextflow.Session
import nextflow.executor.Executor
import nextflow.script.BodyDef
import nextflow.script.ProcessConfig
import spock.lang.Specification
import test.TestHelper
/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class TaskGroupCollectorTest extends Specification {

    def 'should submit tasks as task groups' () {
        given:
        def executor = Mock(Executor)
        def handler = Mock(TaskHandler)
        def taskGroup = [:] as TaskGroup
        def collector = Spy(new TaskGroupCollector(executor, 5)) {
            createTaskGroup(_) >> taskGroup
        }
        and:
        def task = Mock(TaskRun) {
            getConfig() >> Mock(TaskConfig) {
                getAttempt() >> 1
            }
        }

        // collect tasks into task group
        when:
        collector.collect(task)
        collector.collect(task)
        collector.collect(task)
        collector.collect(task)
        then:
        4 * executor.createTaskHandler(task) >> handler
        0 * executor.submit(_)

        // submit task group when it is ready
        when:
        collector.collect(task)
        then:
        1 * executor.createTaskHandler(task) >> handler
        5 * handler.prepareLauncher()
        1 * executor.submit(taskGroup)

        // submit partial task group when closed
        when:
        collector.collect(task)
        collector.close()
        then:
        1 * executor.createTaskHandler(task) >> handler
        1 * handler.prepareLauncher()
        1 * executor.submit(taskGroup)

        // submit tasks directly once closed
        when:
        collector.collect(task)
        then:
        1 * executor.submit(task)
    }

    def 'should submit retried tasks directly' () {
        given:
        def executor = Mock(Executor)
        def collector = Spy(new TaskGroupCollector(executor, 5))
        and:
        def task = Mock(TaskRun) {
            getConfig() >> Mock(TaskConfig) {
                getAttempt() >> 2
            }
        }

        when:
        collector.collect(task)
        then:
        1 * executor.submit(task)
    }

    def 'should create task group' () {
        given:
        def executor = Mock(Executor) {
            getWorkDir() >> TestHelper.createInMemTempDir()
        }
        def collector = Spy(new TaskGroupCollector(executor, 5))
        and:
        def task = Mock(TaskRun) {
            processor >> Mock(TaskProcessor) {
                config >> Mock(ProcessConfig)
                getExecutor() >> executor
                getSession() >> Mock(Session)
                getTaskBody() >> { new BodyDef(null, 'source') }
            }
            getHash() >> HashCode.fromString('0123456789abcdef')
        }
        def handler = Mock(TaskHandler) {
            getTask() >> task
            getWorkDir() >> Paths.get('/work/foo')
            getLaunchCommand() >> ['bash', '-o', 'pipefail', '-c', 'bash /work/foo/.command.run 2>&1 | tee /work/foo/.command.log']
        }

        when:
        def taskGroup = collector.createTaskGroup([handler, handler, handler])
        then:
        taskGroup.config == task.config
        taskGroup.processor == task.processor
        taskGroup.script == '''
            array=( /work/foo /work/foo /work/foo )
            for task_dir in ${array[@]}; do
                export task_dir
                bash -o pipefail -c 'bash ${task_dir}/.command.run 2>&1 | tee ${task_dir}/.command.log' || true
            done
            '''.stripIndent().leftTrim()
        and:
        taskGroup.isContainerEnabled() == false
    }

}