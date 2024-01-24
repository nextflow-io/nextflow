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
import nextflow.executor.TaskArrayAware
import nextflow.script.BodyDef
import nextflow.script.ProcessConfig
import spock.lang.Specification
import test.TestHelper
/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class TaskArrayCollectorTest extends Specification {

    static class DummyExecutor extends Executor implements TaskArrayAware {
        TaskMonitor createTaskMonitor() { null }
        TaskHandler createTaskHandler(TaskRun task) { null }
    }

    def 'should submit tasks as job arrays' () {
        given:
        def executor = Mock(DummyExecutor)
        def handler = Mock(TaskHandler)
        def taskArray = [:] as TaskArray
        def collector = Spy(new TaskArrayCollector(executor, 5)) {
            createTaskArray(_) >> taskArray
        }
        and:
        def task = Mock(TaskRun) {
            getConfig() >> Mock(TaskConfig) {
                getAttempt() >> 1
            }
        }

        // collect tasks into job array
        when:
        collector.collect(task)
        collector.collect(task)
        collector.collect(task)
        collector.collect(task)
        then:
        4 * executor.createTaskHandler(task) >> handler
        0 * executor.submit(_)

        // submit job array when it is ready
        when:
        collector.collect(task)
        then:
        1 * executor.createTaskHandler(task) >> handler
        5 * handler.prepareLauncher()
        1 * executor.submit(taskArray)

        // submit partial job array when closed
        when:
        collector.collect(task)
        collector.close()
        then:
        1 * executor.createTaskHandler(task) >> handler
        1 * handler.prepareLauncher()
        1 * executor.submit(taskArray)

        // submit tasks directly once closed
        when:
        collector.collect(task)
        then:
        1 * executor.submit(task)
    }

    def 'should submit retried tasks directly' () {
        given:
        def executor = Mock(DummyExecutor)
        def collector = Spy(new TaskArrayCollector(executor, 5))
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

    def 'should create task array' () {
        given:
        def executor = Mock(DummyExecutor) {
            getWorkDir() >> TestHelper.createInMemTempDir()
            getArrayIndexName() >> 'ARRAY_JOB_INDEX'
        }
        def collector = Spy(new TaskArrayCollector(executor, 5))
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
        def taskArray = collector.createTaskArray([handler, handler, handler])
        then:
        taskArray.config == task.config
        taskArray.processor == task.processor
        taskArray.script == '''
            array=( /work/foo /work/foo /work/foo )
            export task_dir=${array[ARRAY_JOB_INDEX]}
            bash -o pipefail -c 'bash ${task_dir}/.command.run 2>&1 | tee ${task_dir}/.command.log'
            '''.stripIndent().leftTrim()
        and:
        taskArray.getArraySize() == 3
        taskArray.getContainerConfig().getEnvWhitelist() == [ 'ARRAY_JOB_INDEX' ]
        taskArray.isContainerEnabled() == false
    }

}
