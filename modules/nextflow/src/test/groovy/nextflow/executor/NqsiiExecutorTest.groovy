/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

import nextflow.processor.TaskConfig
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * @author Till Bayer <till.bayer@gmail.com>
 */
class NqsiiExecutorTest extends Specification {

    def testGetCommandLine() {

        given:
        def executor = [:] as NqsiiExecutor
        def task = Mock(TaskRun); task.getName() >> 'hello world'

        expect:
        executor.getSubmitCommandLine(task, Paths.get('/some/path/script.sh') ) == ['qsub', 'script.sh']

    }

    def testHeaders() {

        setup:
        def executor = [:] as NqsiiExecutor

        // mock process
        def proc = Mock(TaskProcessor)

        // task object
        def task = new TaskRun()
        task.processor = proc
        task.workDir = Paths.get('/work/dir')
        task.name = 'task name'

        when:
        task.config = new TaskConfig()
        then:
        executor.getHeaders(task) == '''
                #PBS -N nf-task_name
                #PBS -o /work/dir/.command.log
                #PBS -j o
                #PBS -b 1
                #PBS -l cpunum_job=1
                NXF_CHDIR=/work/dir
                '''
                .stripIndent().leftTrim()


        when:
        task.config = new TaskConfig()
        task.config.queue = 'alpha'
        task.config.time = '1m'
        then:
        executor.getHeaders(task) == '''
                #PBS -N nf-task_name
                #PBS -o /work/dir/.command.log
                #PBS -j o
                #PBS -b 1
                #PBS -q alpha
                #PBS -l cpunum_job=1
                #PBS -l elapstim_req=00:01:00
                NXF_CHDIR=/work/dir
                '''
                .stripIndent().leftTrim()


        when:
        task.config = new TaskConfig()
        task.config.queue = 'alpha'
        task.config.time = '1m'
        task.config.memory = '1m'
        then:
        executor.getHeaders(task) == '''
                #PBS -N nf-task_name
                #PBS -o /work/dir/.command.log
                #PBS -j o
                #PBS -b 1
                #PBS -q alpha
                #PBS -l cpunum_job=1
                #PBS -l elapstim_req=00:01:00
                #PBS -l memsz_job=1mb
                NXF_CHDIR=/work/dir
                '''
                .stripIndent().leftTrim()



        when:
        task.config = new TaskConfig()
        task.config.queue = 'delta'
        task.config.time = '10m'
        task.config.memory = '5m'
        task.config.cpus = 2
        then:
        executor.getHeaders(task) == '''
                #PBS -N nf-task_name
                #PBS -o /work/dir/.command.log
                #PBS -j o
                #PBS -b 1
                #PBS -q delta
                #PBS -l cpunum_job=2
                #PBS -l elapstim_req=00:10:00
                #PBS -l memsz_job=5mb
                NXF_CHDIR=/work/dir
                '''
                .stripIndent().leftTrim()

        when:
        task.config = new TaskConfig()
        task.config.queue = 'delta'
        task.config.time = '1d'
        task.config.memory = '1g'
        task.config.cpus = 8
        then:
        executor.getHeaders(task) == '''
                #PBS -N nf-task_name
                #PBS -o /work/dir/.command.log
                #PBS -j o
                #PBS -b 1
                #PBS -q delta
                #PBS -l cpunum_job=8
                #PBS -l elapstim_req=24:00:00
                #PBS -l memsz_job=1gb
                NXF_CHDIR=/work/dir
                '''
                .stripIndent().leftTrim()

        when:
        task.config = new TaskConfig()
        task.config.queue = 'delta'
        task.config.time = '2d 6h 10m'
        task.config.memory = '2g'
        then:
        executor.getHeaders(task) == '''
                #PBS -N nf-task_name
                #PBS -o /work/dir/.command.log
                #PBS -j o
                #PBS -b 1
                #PBS -q delta
                #PBS -l cpunum_job=1
                #PBS -l elapstim_req=54:10:00
                #PBS -l memsz_job=2gb
                NXF_CHDIR=/work/dir
                '''
                .stripIndent().leftTrim()

    }

    def WorkDirWithBlanks() {

        setup:
        def executor = Spy(NqsiiExecutor)

        // mock process
        def proc = Mock(TaskProcessor)

        // task object
        def task = new TaskRun()
        task.processor = proc
        task.workDir = Paths.get('/work/dir 1')
        task.name = 'task name'

        when:
        task.config = new TaskConfig()
        then:
        executor.getHeaders(task) == '''
                #PBS -N nf-task_name
                #PBS -o "/work/dir\\ 1/.command.log"
                #PBS -j o
                #PBS -b 1
                #PBS -l cpunum_job=1
                NXF_CHDIR=/work/dir\\ 1
                '''
                .stripIndent().leftTrim()

    }

    @Unroll
    def 'should return valid job name given #name'() {
        given:
        def executor = [:] as NqsiiExecutor
        def task = Mock(TaskRun)
        task.getName() >> name

        expect:
        executor.getJobNameFor(task) == expected
        executor.getJobNameFor(task).size() <= 63

        where:
        name        | expected
        'hello'     | 'nf-hello'
        '12 45'     | 'nf-12_45'
        'really-very-very-long-task-name-taking-more-than-63-chars-which-should-be-clipped' | 'nf-really-very-very-long-task-name-taking-more-than-63-chars-wh'

    }

    def testParseJobId() {

        given:
        def executor = [:] as NqsiiExecutor

        expect:
        executor.parseJobId('\nRequest 143993.ace-ssiox submitted to queue: clexpress.\n') == '143993'
        executor.parseJobId('Request 1234.ace-ssiox submitted to queue: queuefoo.') == '1234'
    }


    def testKillTaskCommand() {

        given:
        def executor = [:] as NqsiiExecutor
        expect:
        executor.killTaskCommand('12345') == ['qdel', '12345']

    }

    def testParseQueueStatus() {

        setup:
        def executor = [:] as NqsiiExecutor
        def text =
                """
                RequestID       ReqName  UserName Queue     Pri STT S   Memory      CPU   Elapse R H M Jobs
                --------------- -------- -------- -------- ---- --- - -------- -------- -------- - - - ----
                109441.ace-ssio U2UnCstr userxxx  clbigmem    0 PRR -    8.31G 177109.06   177193 Y Y Y    1
                109447.ace-ssio U7UnCstr userxxx  clbigmem    0 RUN -    6.27G 177152.16   177193 Y Y Y    1
                130945.ace-ssio PRA_28   userxxx  clbigmem    0 STG -   36.57G 1393093.29    94243 Y Y Y    1
                131687.ace-ssio Decon2   userxxx  clbigmem    0 CHK -   16.69G 29977.38    30013 Y Y Y    1
                131689.ace-ssio Decon1   userxxx  clbigmem    0 EXT -    0.00B     0.00        0 Y Y Y    1
                141025.ace-ssio Chlorell userxxx  clbigmem    0 POR -    0.00B     0.00        0 Y Y Y    1
                142867.ace-ssio butter   userxxx  clexpres    0 ARI -    2.81G 56568.09     3763 Y Y Y    1
                143839.ace-ssio A11148.g userxxx  clexpres    0 FWD -    0.00B     0.00        0 Y Y Y    1
                142611.ace-ssio RUN_Wb6  userxxx  smallque    0 MIG -    5.42G   427.44     9261 Y Y Y    2
                142616.ace-ssio RUN_CT37 userxxx  smallque    0 QUE -    5.42G  1102.06     9230 Y Y Y    2
                142699.ace-ssio RUN_CT44 userxxx  smallque    0 WAT -    5.06G    73.94     8781 Y Y Y    2
                142714.ace-ssio RUN_CT43 userxxx  smallque    0 GQD -   79.06M    10.94     8661 Y Y Y    2
                143311.ace-ssio RUN_CT39 userxxx  smallque    0 SUS -   83.02M     7.25     5181 Y Y Y    2
                143527.ace-ssio RUN_S07  userxxx  smallque    0 HLD -    6.24G   880.60     4371 Y Y Y    2
                143934.ace-ssio RUN_P28  userxxx  smallque    0 HOL -   36.94G 66930.90     1911 Y Y Y    9

                """.stripIndent().trim()

        when:
        def result = executor.parseQueueStatus(text)
        then:
        result.size() == 15
        result['109441'] == AbstractGridExecutor.QueueStatus.RUNNING
        result['109447'] == AbstractGridExecutor.QueueStatus.RUNNING
        result['130945'] == AbstractGridExecutor.QueueStatus.RUNNING
        result['131687'] == AbstractGridExecutor.QueueStatus.RUNNING
        result['131689'] == AbstractGridExecutor.QueueStatus.RUNNING
        result['141025'] == AbstractGridExecutor.QueueStatus.RUNNING
        result['142867'] == AbstractGridExecutor.QueueStatus.PENDING
        result['143839'] == AbstractGridExecutor.QueueStatus.PENDING
        result['142611'] == AbstractGridExecutor.QueueStatus.PENDING
        result['142616'] == AbstractGridExecutor.QueueStatus.PENDING
        result['142699'] == AbstractGridExecutor.QueueStatus.PENDING
        result['142714'] == AbstractGridExecutor.QueueStatus.PENDING
        result['143311'] == AbstractGridExecutor.QueueStatus.HOLD
        result['143527'] == AbstractGridExecutor.QueueStatus.HOLD
        result['143934'] == AbstractGridExecutor.QueueStatus.HOLD

    }

    def 'should return qstat command line' () {
        given:
        def executor = [:] as NqsiiExecutor

        expect:
        executor.queueStatusCommand(null) == ['qstat']
        executor.queueStatusCommand('queuexxx') == ['qstat', '-q', 'queuexxx']
        executor.queueStatusCommand('xxx').each { assert it instanceof String }
    }

}
