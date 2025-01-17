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
            1 FAILED
            2 FINISHED
            3 RUNNING
            '''.stripIndent().leftTrim()
        when:
        def result = exec.parseQueueStatus(STATUS)
        then:
        result.size() == 3
        and:
        result.get('1') == AbstractGridExecutor.QueueStatus.ERROR
        result.get('2') == AbstractGridExecutor.QueueStatus.DONE
        result.get('3') == AbstractGridExecutor.QueueStatus.RUNNING
    }

    def 'should parse job id' () {
        given:
        def exec = new HyperQueueExecutor()
        and:
        def JOB_OK = '''
            3
            '''.stripIndent().leftTrim()

        def JOB_FAIL = '''
            error: submit failed this and that
            '''.stripIndent().leftTrim()

        
        when:
        def result = exec.parseJobId(JOB_OK)
        then:
        result == '3'

        when:
        exec.parseJobId(JOB_FAIL)
        then:
        def err = thrown(IllegalArgumentException)
        err.message == '''
            HyperQueue submit failed or response is invalid:
            error: submit failed this and that


            '''.stripIndent().leftTrim()
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
        result == ['hq', '--output-mode=quiet', 'submit', '--directives=file', 'script.run']
    }

    def 'should get kill command' () {
        when:
        def executor = Spy(HyperQueueExecutor)
        then:
        executor.killTaskCommand('12345').join(' ') == 'hq job cancel 12345'
        executor.killTaskCommand(['12345','12']).join(' ') == 'hq job cancel 12345,12'
    }

    def 'should validate headers' () {

        setup:
        def executor = new HyperQueueExecutor()
        def proc = Mock(TaskProcessor)
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
        task.config.accelerator = [request: 1, limit: 2]
        then:
        executor.getHeaders(task) == '''
            #HQ --name nf-task-1
            #HQ --log /work/dir/.command.log
            #HQ --cwd /work/dir
            #HQ --resource mem=8192
            #HQ --cpus 4
            #HQ --time-limit 60sec
            #HQ --resource gpus=2
            '''
            .stripIndent().leftTrim()
    }
}
