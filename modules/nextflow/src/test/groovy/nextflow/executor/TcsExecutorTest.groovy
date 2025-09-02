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
 */

package nextflow.executor

import java.nio.file.Paths

import nextflow.Session
import nextflow.processor.TaskArrayRun
import nextflow.processor.TaskConfig
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Satoshi Ohshima <ohshima@cc.kyushu-u.ac.jp>
 */

class TcsExecutorTest extends Specification {

    def testParseJob() {
        given:
        def exec = [:] as TcsExecutor

        expect:
        exec.parseJobId('[INFO] PJM 0000 pjsub Job 123456 submitted.') == '123456'
        exec.parseJobId('[INFO] PJM 0000 pjsub Job 630114 submitted.') == '630114'
    }

    def testKill() {
        given:
        def exec = [:] as TcsExecutor

        expect:
        exec.killTaskCommand(123) == ['pjdel','123']
    }

    @Unroll
    def testGetCommandLine() {
        given:
        def session = Mock(Session) {getConfig()>>[:]}
        def exec = Spy(TcsExecutor) { getSession()>>session }

        when:
        def result = exec.getSubmitCommandLine(Mock(TaskRun), Paths.get(PATH))
        then:
        exec.pipeLauncherScript() >> PIPE
        result == EXPECTED

        where:
        PATH                    | PIPE      | EXPECTED
        '/some/path/job.sh'     | false     | ['pjsub', 'job.sh']
        '/some/path/job.sh'     | true      | ['pjsub']
    }

    def 'test job script headers' () {
        setup:
        // TCS executor
        def executor = [:] as TcsExecutor
        executor.session = Mock(Session)

        // mock process
        def proc = Mock(TaskProcessor)

        // task object
        def task = new TaskRun()
        task.processor = proc
        task.workDir = Paths.get('/work/path')
        task.name = 'the task name'

        // 1min job
        when:
		task.index = 21
        task.config = new TaskConfig()
		task.config.time = '00:10:00'
        then:
        executor.getHeaders(task) == '''
                #PJM -N nf-the_task_name
                #PJM -L elapse=00:10:00
                #PJM -o /work/path/.command.log
                #PJM -j
                #PJM -S
                '''
                .stripIndent().leftTrim()

        // 1hour job
        when:
		task.index = 22
        task.config = new TaskConfig()
		task.config.time = '01:00:00'
        task.name = '1hour job'
        then:
        executor.getHeaders(task) == '''
                #PJM -N nf-1hour_job
                #PJM -L elapse=01:00:00
                #PJM -o /work/path/.command.log
                #PJM -j
                #PJM -S
                '''
                .stripIndent().leftTrim()

        // job array (bulk job)
        when: 'with job array (bulk job)'
		task.config = new TaskConfig()
		task.workDir = Paths.get('/work/path')
		task.name = 'array job'
        def taskArray = Mock(TaskArrayRun) {
            config >> new TaskConfig()
            name >> task.name
            getArraySize() >> 5
            workDir >> task.workDir
        }
		taskArray.config.time = '01:00:00'

        then:
        executor.getHeaders(taskArray) == '''
                #PJM -N nf-array_job
                #PJM -L elapse=01:00:00
                #PJM -j
                #PJM -S
                '''
                .stripIndent().leftTrim()
    }

    def testWorkDirWithBlanks() {
        setup:
        // LSF executor
        def executor = Spy(TcsExecutor)
        executor.session = Mock(Session)

        // mock process
        def proc = Mock(TaskProcessor)

        // task object
        def task = new TaskRun()
        task.processor = proc
        task.workDir = Paths.get('/work/some data/path')
        task.name = 'the task name'

        when:
        task.index = 21
        task.config = new TaskConfig()
		task.config.time = '00:10:00'
        then:
        executor.getHeaders(task) == '''
                #PJM -N nf-the_task_name
                #PJM -L elapse=00:10:00
                #PJM -o "/work/some\\ data/path/.command.log"
                #PJM -j
                #PJM -S
                '''
                .stripIndent().leftTrim()
    }

    def testQstatCommand() {
        given:
        def executor = [:] as TcsExecutor
        def text =
                """
                100001    job1      NM ACC
                100002    job2      NM QUE
                100003    job3      NM RNA
                100004    job4      NM RUN
                100005    job5      NM RNO
                100006    job6      NM EXT
                100007    job7      NM CCL
                100008    job8      NM HLD
                100009    job9      NM ERR
                """.stripIndent().trim()


        when:
        def result = executor.parseQueueStatus(text)
        then:
        result.size() == 9
        result['100001'] == AbstractGridExecutor.QueueStatus.PENDING
        result['100002'] == AbstractGridExecutor.QueueStatus.PENDING
        result['100003'] == AbstractGridExecutor.QueueStatus.RUNNING
        result['100004'] == AbstractGridExecutor.QueueStatus.RUNNING
        result['100005'] == AbstractGridExecutor.QueueStatus.RUNNING
        result['100006'] == AbstractGridExecutor.QueueStatus.DONE
        result['100007'] == AbstractGridExecutor.QueueStatus.DONE
        result['100008'] == AbstractGridExecutor.QueueStatus.HOLD
        result['100009'] == AbstractGridExecutor.QueueStatus.ERROR
    }

    def testQueueStatusCommand() {
        when:
        def exec = [:] as TcsExecutor
        then:
        exec.queueStatusCommand(null) == ['pjstat', '-E']
    }

    def 'should get array index name and start' () {
        given:
        def executor = Spy(TcsExecutor)
        expect:
        executor.getArrayIndexName() == 'PJM_BULKNUM'
        executor.getArrayIndexStart() == 0
    }

    @Unroll
    def 'should get array task id' () {
        given:
        def executor = Spy(TcsExecutor)
        expect:
        executor.getArrayTaskId(JOB_ID, TASK_INDEX) == EXPECTED

        where:
        JOB_ID      | TASK_INDEX    | EXPECTED
        '1234'      | 1             | '1234[1]'
        '123456'    | 2             | '123456[2]'
    }

}

