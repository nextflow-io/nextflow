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
import nextflow.script.BaseScript
import nextflow.script.BodyDef
import nextflow.script.ProcessConfig
import spock.lang.Specification
import test.TestHelper
/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class TaskBatchCollectorTest extends Specification {

    def 'should submit tasks as task batches' () {
        given:
        def executor = Mock(Executor)
        def handler = Mock(TaskHandler)
        def taskBatch = [:] as TaskBatchRun
        def collector = Spy(new TaskBatchCollector(null, executor, 5, false)) {
            createTaskBatch(_) >> taskBatch
        }
        and:
        def task = Mock(TaskRun) {
            getConfig() >> Mock(TaskConfig) {
                getAttempt() >> 1
            }
        }

        // collect tasks into task batch
        when:
        collector.collect(task)
        collector.collect(task)
        collector.collect(task)
        collector.collect(task)
        then:
        0 * executor.submit(_)

        // submit task batch when it is ready
        when:
        collector.collect(task)
        then:
        1 * executor.submit(taskBatch)

        // submit partial task batch when closed
        when:
        collector.collect(task)
        collector.collect(task)
        collector.close()
        then:
        1 * executor.submit(taskBatch)

        // submit tasks directly once closed
        when:
        collector.collect(task)
        then:
        1 * executor.submit(task)
    }

    def 'should submit retried tasks directly' () {
        given:
        def executor = Mock(Executor)
        def collector = Spy(new TaskBatchCollector(null, executor, 5, false))
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

    def 'should create task batch' () {
        given:
        def executor = Mock(Executor) {
            getWorkDir() >> TestHelper.createInMemTempDir()
        }
        def config = Spy(ProcessConfig, constructorArgs: [Mock(BaseScript), 'PROC']) {
            createTaskConfig() >> Mock(TaskConfig)
            get('cpus') >> 4
            get('tag') >> 'foo'
        }
        def proc = Mock(TaskProcessor) {
            getConfig() >> config
            getExecutor() >> executor
            getName() >> 'PROC'
            getSession() >> Mock(Session)
            isSingleton() >> false
            getTaskBody() >> { new BodyDef(null, 'source') }
        }
        def collector = Spy(new TaskBatchCollector(proc, executor, 5, false))
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
        def taskBatch = collector.createTaskBatch([task, task, task])
        then:
        3 * executor.createTaskHandler(task) >> handler
        3 * handler.prepareLauncher()
        1 * collector.createBatchTaskScript([handler, handler, handler]) >> 'the-task-batch-script'
        and:
        taskBatch.name == 'PROC (1)'
        taskBatch.config.cpus == 4
        taskBatch.config.tag == null
        taskBatch.processor == proc
        taskBatch.script == 'the-task-batch-script'
        and:
        taskBatch.children.size() == 3
        taskBatch.isContainerEnabled() == false
    }

    def 'should get batch launch script' () {
        given:
        def processor = Mock(TaskProcessor)
        def executor = Spy(Executor) {
            isWorkDirDefaultFS() >> true
        }
        def collector = Spy(new TaskBatchCollector(processor, executor, 10, false))
        and:
        def h1 = Mock(TaskHandler)
        def h2 = Mock(TaskHandler)
        def h3 = Mock(TaskHandler)

        when:
        def result = collector.createBatchTaskScript([h1,h2,h3])
        then:
        executor.getChildWorkDir(h1) >> '/work/dir/1'
        executor.getChildWorkDir(h2) >> '/work/dir/2'
        executor.getChildWorkDir(h3) >> '/work/dir/3'
        then:
        result == '''\
            array=( /work/dir/1 /work/dir/2 /work/dir/3 )
            for nxf_batch_task_dir in ${array[@]}; do
                export nxf_batch_task_dir
                bash $nxf_batch_task_dir/.command.run 2>&1 > $nxf_batch_task_dir/.command.log || true
            done

            '''.stripIndent(true)
    }

    def 'should execute tasks in parallel if specified' () {
        given:
        def processor = Mock(TaskProcessor)
        def executor = Spy(Executor) {
            isWorkDirDefaultFS() >> true
        }
        def collector = Spy(new TaskBatchCollector(processor, executor, 10, true))
        and:
        def h1 = Mock(TaskHandler)
        def h2 = Mock(TaskHandler)
        def h3 = Mock(TaskHandler)

        when:
        def result = collector.createBatchTaskScript([h1,h2,h3])
        then:
        executor.getChildWorkDir(h1) >> '/work/dir/1'
        executor.getChildWorkDir(h2) >> '/work/dir/2'
        executor.getChildWorkDir(h3) >> '/work/dir/3'
        then:
        result == '''\
            array=( /work/dir/1 /work/dir/2 /work/dir/3 )
            for nxf_batch_task_dir in ${array[@]}; do
                export nxf_batch_task_dir
                bash $nxf_batch_task_dir/.command.run 2>&1 > $nxf_batch_task_dir/.command.log || true &
            done
            wait
            '''.stripIndent(true)
    }

}