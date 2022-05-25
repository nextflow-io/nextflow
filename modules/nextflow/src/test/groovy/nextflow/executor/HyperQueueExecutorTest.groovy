/*
 * Copyright 2020-2022, Seqera Labs
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

package nextflow.executor

import java.nio.file.Paths

import nextflow.processor.TaskConfig
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import spock.lang.Specification

/**
 * Test HyperQueue executor
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * @author Henrik Nortamo <henrik.nortamo@csc.fi>
 */
class HyperQueueExecutorTest extends Specification {

    def 'should parse status' () {
        given:
        def exec = new HyperQueueExecutor()
        and:
        def STATUS = '''
        [
          {
            "id": 1,
            "name": "tash.sk",
            "task_count": 1,
            "task_stats": {
              "canceled": 0,
              "failed": 1,
              "finished": 0,
              "running": 0,
              "waiting": 0
            }
          },
          {
            "id": 2,
            "name": "task.sh",
            "task_count": 1,
            "task_stats": {
              "canceled": 0,
              "failed": 0,
              "finished": 1,
              "running": 0,
              "waiting": 0
            }
          },
          {
            "id": 3,
            "name": "task.sh",
            "task_count": 1,
            "task_stats": {
              "canceled": 0,
              "failed": 0,
              "finished": 0,
              "running": 1,
              "waiting": 0
            }
          }
        ]
        '''
        when:
        def result = exec.parseQueueStatus(STATUS)
        then:
        result.size() == 3
        and:
        result.get('1') == AbstractGridExecutor.QueueStatus.ERROR
        result.get('2') == AbstractGridExecutor.QueueStatus.DONE
        result.get('3') == AbstractGridExecutor.QueueStatus.RUNNING
    }

    def 'should decode status' () {
        given:
        def exec = new HyperQueueExecutor()

        expect:
        exec.parseJobStatus([task_stats: [failed: 1]]) == AbstractGridExecutor.QueueStatus.ERROR
        exec.parseJobStatus([task_stats: [canceled: 1]]) == AbstractGridExecutor.QueueStatus.DONE
        exec.parseJobStatus([task_stats: [finished: 1]]) == AbstractGridExecutor.QueueStatus.DONE
        exec.parseJobStatus([task_stats: [running: 1]]) == AbstractGridExecutor.QueueStatus.RUNNING
        exec.parseJobStatus([task_stats: [waiting: 1]]) == AbstractGridExecutor.QueueStatus.HOLD
        exec.parseJobStatus([:]) == AbstractGridExecutor.QueueStatus.UNKNOWN
    }

    def 'should parse job id' () {
        given:
        def exec = new HyperQueueExecutor()
        and:
        def JOB_OK = '''
        {
            "id": 3
        }
        '''

        def JOB_FAIL = '''
        {
            "error": "submit failed this and that"
        }
        '''

        
        when:
        def result = exec.parseJobId(JOB_OK)
        then:
        result == '3'

        when:
        exec.parseJobId(JOB_FAIL)
        then:
        def err = thrown(IllegalArgumentException)
        err.message == 'HyperQueue submit error: submit failed this and that'

        when:
        exec.parseJobId('not ok')
        then:
        err = thrown(IllegalArgumentException)
        err.message == '''\
            Invalid HyperQueue submit JSON response:
            not ok
                    
            '''.stripIndent()
    }

    def 'should create submit command' () {
        given:
        def exec = new HyperQueueExecutor()
        and:
        def task = Mock(TaskRun)
        def path = Paths.get('/some/work/script.run')

        when:
        def result = exec.getSubmitCommandLine(task, path)
        then:
        result == ['hq', '--output-mode=json','submit', '--directives=file', 'script.run']
    }

    def 'should get kill command'() {
        when:
        // executor stub object
        def executor = Spy(HyperQueueExecutor)
        then:
        executor.killTaskCommand('12345').join(' ') == 'hq job cancel 12345'
        executor.killTaskCommand(['12345','12']).join(' ') == 'hq job cancel 12345,12'

    }

    def "should validate headers"() {

        setup:
        def executor = new HyperQueueExecutor()

        // mock process
        def proc = Mock(TaskProcessor)

        // task object
        def task = new TaskRun()
        task.processor = proc
        task.workDir = Paths.get('/work/dir')
        task.name = 'task-1'

        when:
        task.config = new TaskConfig()
        then:
        executor.getHeaders(task) == '''
                #HQ --name nf-task-1
                #HQ --log /work/dir/.command.log
                #HQ --cwd /work/dir
                '''
                .stripIndent().leftTrim()


        when:
        task.config = new TaskConfig()
        task.config.cpus = 4
        task.config.time = '1m'
        task.config.memory = '8GB'
        task.config.accelerator = [request: 1, limit: 2, type: 'nvida']
        then:
        executor.getHeaders(task) == '''
                #HQ --name nf-task-1
                #HQ --log /work/dir/.command.log
                #HQ --cwd /work/dir
                #HQ --resource mem=8589934592
                #HQ --cpus 4
                #HQ --time-limit 60sec
                #HQ --resource gpus=2
                '''
                .stripIndent().leftTrim()
    }
}
