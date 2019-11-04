/*
 * Copyright 2019, Hospices Civils de Lyon (HCL)
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
/**
 *
 * @author Maxime Vallée <maxime.vallee@chu-lyon.fr>
 */
class OarExecutorTest extends Specification {

    def testParseJob() {

        given:
        def exec = [:] as OarExecutor

        expect:
        exec.parseJobId('OAR_JOB_ID=10') == '10'
        exec.parseJobId('OAR_JOB_ID=20') == '20'
        exec.parseJobId('OAR_JOB_ID=30') == '30'
        exec.parseJobId('OAR_JOB_ID=40\n') == '40'
        exec.parseJobId('\nOAR_JOB_ID=50') == '50'

        when:
        exec.parseJobId('Something else 10')
        then:
        thrown(IllegalStateException)

    }

    def testKill() {

        given:
        def exec = [:] as OarExecutor
        expect:
        exec.killTaskCommand(123) == ['oardel','123']

    }

    def testGetCommandLine() {

        when:
        def exec = [:] as OarExecutor
        then:
        exec.getSubmitCommandLine(Mock(TaskRun), Paths.get('/some/path/job.sh')) == ['oarsub', '-S', './job.sh']
    }

    def testGetHeaders() {

        setup:
        // OAR executor
        def executor = [:] as OarExecutor

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
        then:
        executor.getHeaders(task) == '''
                #OAR -d /work/path
                #OAR -n nf-the_task_name
                #OAR -O /work/path/.command.out
                #OAR -E /work/path/.command.err
                '''
                .stripIndent().leftTrim()

        when:
        task.config = new TaskConfig()
        task.config.queue = 'delta'
        then:
        executor.getHeaders(task) == '''
                #OAR -d /work/path
                #OAR -n nf-the_task_name
                #OAR -O /work/path/.command.out
                #OAR -E /work/path/.command.err
                #OAR -q delta
                '''
                .stripIndent().leftTrim()

        when:
        task.config = new TaskConfig()
        task.config.time = '1m'
        then:
        executor.getHeaders(task) == '''
                #OAR -d /work/path
                #OAR -n nf-the_task_name
                #OAR -O /work/path/.command.out
                #OAR -E /work/path/.command.err
                #OAR -l walltime=00:01:00
                '''
                .stripIndent().leftTrim()

        when:
        task.config = new TaskConfig()
        task.config.time = '1h'
        task.config.memory = '1023 M'
        task.config.clusterOptions = '-a 1'
        then:
        executor.getHeaders(task) == '''
                #OAR -d /work/path
                #OAR -n nf-the_task_name
                #OAR -O /work/path/.command.out
                #OAR -E /work/path/.command.err
                #OAR -p "memnode=1"
                #OAR -l walltime=01:00:00
                #OAR -a 1
                '''
                .stripIndent().leftTrim()

        when:
        task.config = new TaskConfig()
        task.config.time = '1h'
        task.config.memory = '2049 M'
        task.config.clusterOptions = '-a 1'
        then:
        executor.getHeaders(task) == '''
                #OAR -d /work/path
                #OAR -n nf-the_task_name
                #OAR -O /work/path/.command.out
                #OAR -E /work/path/.command.err
                #OAR -p "memnode=2"
                #OAR -l walltime=01:00:00
                #OAR -a 1
                '''
                .stripIndent().leftTrim()

        when:
        task.config = new TaskConfig()
        task.config.cpus = 2
        task.config.time = '2h'
        task.config.memory = '4 G'
        task.config.clusterOptions = '-b 2'
        then:
        executor.getHeaders(task) == '''
                #OAR -d /work/path
                #OAR -n nf-the_task_name
                #OAR -O /work/path/.command.out
                #OAR -E /work/path/.command.err
                #OAR -p "memnode=4"
                #OAR -l /nodes=1/core=2,walltime=02:00:00
                #OAR -b 2
                '''
                .stripIndent().leftTrim()

        when:
        task.config = new TaskConfig()
        task.config.cpus = 8
        task.config.time = '2d 3h'
        task.config.memory = '6 G'
        task.config.clusterOptions = '-x 3'
        then:
        executor.getHeaders(task) == '''
                #OAR -d /work/path
                #OAR -n nf-the_task_name
                #OAR -O /work/path/.command.out
                #OAR -E /work/path/.command.err
                #OAR -p "memnode=6"
                #OAR -l /nodes=1/core=8,walltime=51:00:00
                #OAR -x 3
                '''
                .stripIndent().leftTrim()
    }

    def testWorkDirWithBlanks() {

        setup:
        // OAR executor
        def executor = Spy(OarExecutor)

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
        then:
        executor.getHeaders(task) == '''
                #OAR -d "/work/some\\ data/path"
                #OAR -n nf-the_task_name
                #OAR -O "/work/some\\ data/path/.command.out"
                #OAR -E "/work/some\\ data/path/.command.err"
                '''
                .stripIndent().leftTrim()

    }


    def testQstatCommand() {

        setup:
        def executor = [:] as OarExecutor
        def text =
                """
                Job_Id: 4930950
                    state = toLaunch
                Job_Id: 4930951
                    state = Launching
                Job_Id: 4930952
                    state = Running
                Job_Id: 4930953
                    state = Finishing
                Job_Id: 4930954
                    state = Waiting
                Job_Id: 4930955
                    state = toAckReservation
                Job_Id: 4930956
                    state = Hold
                Job_Id: 4930957
                    state = Suspended
                Job_Id: 4930958
                    state = Error
                Job_Id: 4930959
                    state = toError
                Job_Id: 4930960
                    state = Terminated
                """.stripIndent().trim()


        when:
        def result = executor.parseQueueStatus(text)
        then:
        result.size() == 11
        result['4930950'] == AbstractGridExecutor.QueueStatus.RUNNING
        result['4930951'] == AbstractGridExecutor.QueueStatus.RUNNING
        result['4930952'] == AbstractGridExecutor.QueueStatus.RUNNING
        result['4930953'] == AbstractGridExecutor.QueueStatus.RUNNING
        result['4930954'] == AbstractGridExecutor.QueueStatus.PENDING
        result['4930955'] == AbstractGridExecutor.QueueStatus.PENDING
        result['4930956'] == AbstractGridExecutor.QueueStatus.HOLD
        result['4930957'] == AbstractGridExecutor.QueueStatus.HOLD
        result['4930958'] == AbstractGridExecutor.QueueStatus.ERROR
        result['4930959'] == AbstractGridExecutor.QueueStatus.ERROR
        result['4930960'] == AbstractGridExecutor.QueueStatus.DONE

    }

    def testQueueStatusCommand() {
        when:
        def executor = [:] as OarExecutor
        then:
        executor.queueStatusCommand(null) == ['sh','-c', "oarstat -f | egrep '(Job_Id:|state =)' ".toString()]
        executor.queueStatusCommand('xxx') == ['sh','-c', "oarstat -f xxx | egrep '(Job_Id:|state =)' ".toString()]

    }

}
