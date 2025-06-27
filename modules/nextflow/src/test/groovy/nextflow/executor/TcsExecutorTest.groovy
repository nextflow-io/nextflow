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
 * @author Satoshi Ohshima <xxx@xxx.xxxx>
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

        when:
		task.index = 21
        task.config = new TaskConfig()
		task.config.time = '00:10:00'
        then:
        executor.getHeaders(task) == '''
		        #PJM -N nf-the_task_name
				#PJM -j
				#PJM -S
				#PJM -L elapse=00:10:00
                '''
                .stripIndent().leftTrim()
    }
/*
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
				#PJM -j
				#PJM -S
				#PJM -L elapse=00:10:00
                #PJM -o "/work/some\\ data/path/.command.log"
                '''
                .stripIndent().leftTrim()

    }
*/

    def testQstatCommand() {

        setup:
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
        def usr = System.getProperty('user.name')
        def exec = [:] as TcsExecutor
        then:
        usr
        exec.queueStatusCommand(null) == ['pjstat2', '-E']
        //exec.queueStatusCommand('long') == ['pjstat2','-E']
    }

    def 'should get array (bulk) index name and start' () {
        given:
        def executor = Spy(TcsExecutor)
        expect:
        executor.getArrayIndexName() == 'PJM_BULKNUM'
        executor.getArrayIndexStart() == 0
    }

    @Unroll
    def 'should get array (bulk) task id' () {
        given:
        def executor = Spy(TcsExecutor)
        expect:
        executor.getArrayTaskId(JOB_ID, TASK_INDEX) == EXPECTED

        where:

        JOB_ID      | TASK_INDEX    | EXPECTED
        '1234'      | 1             | '1234[1]'
        '123456'    | 2             | '123456[2]'
/*
        JOB_ID      JOB_NAME  MD  STATUS  USER        RSCGROUP     START_DATE             ELAPSE    NODE  CORE  GPU  POINT
        3209914     bulk1.sh  BU  RUN     ku40000105  a-batch-low  -                      -            -     1    -      -
        3209914[1]  bulk1.sh  BU  RUN     ku40000105  a-batch-low  2025/06/13 17:53:43    00:00:01     -     1    -      1
*/
    }

/* There is no project option in TCS. (on Genkai and Flow) */
/*
    @Unroll
    def 'should set tcs account' () {
        given:
        // task
        def task = new TaskRun()
        task.workDir = Paths.get('/work/dir')
        task.processor = Mock(TaskProcessor)
        task.processor.getSession() >> Mock(Session)
        task.config = Mock(TaskConfig)
        and:
        def executor = Spy(TcsExecutor)
        executor.getJobNameFor(_) >> 'foo'
        executor.getName() >> 'tcs'
        executor.getSession() >> Mock(Session) { getExecConfigProp('tcs', 'account',null)>>ACCOUNT }

        when:
        def result = executor.getDirectives(task, [])
        then:
        result == EXPECTED

        where:
        ACCOUNT             | EXPECTED
        null                | ['-N', 'foo', '-o', '/work/dir/.command.log', '-j']
        'project-123'       | ['-N', 'foo', '-o', '/work/dir/.command.log', '-j']
    }
	*/
}
