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
import nextflow.executor.TaskArrayExecutor
import nextflow.script.BaseScript
import nextflow.script.BodyDef
import nextflow.script.ProcessConfig
import spock.lang.Specification
import test.TestHelper
/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class TaskArrayCollectorTest extends Specification {

    static class DummyExecutor extends Executor implements TaskArrayExecutor {
        TaskMonitor createTaskMonitor() { null }
        TaskHandler createTaskHandler(TaskRun task) { null }

        String getArrayIndexName() { null }
        int getArrayIndexStart() { 0 }
        String getArrayTaskId(String jobId, int index) { null }
    }

    def 'should submit tasks as job arrays' () {
        given:
        def executor = Mock(DummyExecutor)
        def handler = Mock(TaskHandler)
        def taskArray = [:] as TaskArrayRun
        def collector = Spy(new TaskArrayCollector(null, executor, 5)) {
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
        0 * executor.submit(_)

        // submit job array when it is ready
        when:
        collector.collect(task)
        then:
        1 * executor.submit(taskArray)

        // submit partial job array when closed
        when:
        collector.collect(task)
        collector.collect(task)
        collector.close()
        then:
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
        def collector = Spy(new TaskArrayCollector(null, executor, 5))
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
        def exec = Mock(DummyExecutor) {
            getWorkDir() >> TestHelper.createInMemTempDir()
            getArrayIndexName() >> 'ARRAY_JOB_INDEX'
        }
        def config = Spy(ProcessConfig, constructorArgs: [Mock(BaseScript), 'PROC']) {
            createTaskConfig() >> Mock(TaskConfig)
            get('cpus') >> 4
            get('tag') >> 'foo'
        }
        def proc = Mock(TaskProcessor) {
            getConfig() >> config
            getExecutor() >> exec
            getName() >> 'PROC'
            getSession() >> Mock(Session)
            isSingleton() >> false
            getTaskBody() >> { new BodyDef(null, 'source') }
        }
        def collector = Spy(new TaskArrayCollector(proc, exec, 5))
        and:
        def task = Mock(TaskRun) {
            index >> 1
            getHash() >> HashCode.fromString('0123456789abcdef')
            getWorkDir() >> Paths.get('/work/foo')
        }
        def handler = Mock(TaskHandler) {
            getTask() >> task
        }

        when:
        def taskArray = collector.createTaskArray([task, task, task])
        then:
        3 * exec.createTaskHandler(task) >> handler
        3 * handler.prepareLauncher()
        1 * collector.createArrayTaskScript([handler, handler, handler]) >> 'the-task-array-script'
        and:
        taskArray.name == 'PROC (1)'
        taskArray.config.cpus == 4
        taskArray.config.tag == null
        taskArray.processor == proc
        taskArray.script == 'the-task-array-script'
        and:
        taskArray.getArraySize() == 3
        taskArray.getContainerConfig().getEnvWhitelist() == [ 'ARRAY_JOB_INDEX' ]
        taskArray.isContainerEnabled() == false
    }

    def 'should get array ref' () {
        given:
        def processor = Mock(TaskProcessor)
        def executor = Mock(DummyExecutor) {
            getArrayIndexName() >> NAME
            getArrayIndexStart() >> START
        }
        def collector = Spy(new TaskArrayCollector(processor, executor, 10))
        expect:
        collector.getArrayIndexRef() == EXPECTED

        where:
        NAME        | START     | EXPECTED
        'INDEX'     | 0         | '${array[INDEX]}'
        'INDEX'     | 1         | '${array[INDEX - 1]}'
        'OTHER'     | 99        | '${array[OTHER - 99]}'
    }

    def 'should get array launch script' () {
        given:
        def processor = Mock(TaskProcessor)
        def executor = Spy(DummyExecutor) { isWorkDirDefaultFS()>>true }
        def collector = Spy(new TaskArrayCollector(processor, executor, 10))
        and:
        def h1 = Mock(TaskHandler)
        def h2 = Mock(TaskHandler)
        def h3 = Mock(TaskHandler)

        when:
        def result = collector.createArrayTaskScript([h1,h2,h3])
        then:
        executor.getArrayWorkDir(h1) >> '/work/dir/1'
        executor.getArrayWorkDir(h2) >> '/work/dir/2'
        executor.getArrayWorkDir(h3) >> '/work/dir/3'
        and:
        collector.getArrayIndexRef() >> '$array[INDEX]'
        then:
        result == '''\
                array=( /work/dir/1 /work/dir/2 /work/dir/3 )
                export nxf_array_task_dir=$array[INDEX]
                bash $nxf_array_task_dir/.command.run 2>&1 > $nxf_array_task_dir/.command.log
                '''.stripIndent(true)
    }
}
