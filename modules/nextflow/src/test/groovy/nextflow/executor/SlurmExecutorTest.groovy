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
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class SlurmExecutorTest extends Specification {

    def testParseJob() {

        given:
        def exec = [:] as SlurmExecutor

        expect:
        exec.parseJobId('Submitted batch job 10') == '10'
        exec.parseJobId('Submitted batch job 20') == '20'
        exec.parseJobId('30') == '30'
        exec.parseJobId('40\n') == '40'
        exec.parseJobId('\n50') == '50'
        exec.parseJobId('Submitted batch job 630114 on cluster mpp2') == '630114'

        when:
        exec.parseJobId('Something else 10')
        then:
        thrown(IllegalStateException)

    }

    def testKill() {

        given:
        def exec = [:] as SlurmExecutor
        expect:
        exec.killTaskCommand(123) == ['scancel','123']

    }

    @Unroll
    def testGetCommandLine() {
        given:
        def session = Mock(Session) {getConfig()>>[:]}
        def exec = Spy(SlurmExecutor) { getSession()>>session }

        when:
        def result = exec.getSubmitCommandLine(Mock(TaskRun), Paths.get(PATH))
        then:
        exec.pipeLauncherScript() >> PIPE
        result == EXPECTED

        where:
        PATH                    | PIPE      | EXPECTED
        '/some/path/job.sh'     | false     | ['sbatch', 'job.sh']
        '/some/path/job.sh'     | true      | ['sbatch']
    }

    def 'test job script headers' () {

        setup:
        // SLURM executor
        def executor = [:] as SlurmExecutor
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
        then:
        executor.getHeaders(task) == '''
                #SBATCH -J nf-the_task_name
                #SBATCH -o /work/path/.command.log
                #SBATCH --no-requeue
                #SBATCH --signal B:USR2@30
                '''
                .stripIndent().leftTrim()

        when:
        task.config = new TaskConfig()
        task.config.queue = 'delta'
        then:
        executor.getHeaders(task) == '''
                #SBATCH -J nf-the_task_name
                #SBATCH -o /work/path/.command.log
                #SBATCH --no-requeue
                #SBATCH --signal B:USR2@30
                #SBATCH -p delta
                '''
                .stripIndent().leftTrim()

        when:
        task.config = new TaskConfig()
        task.config.time = '1m'
        then:
        executor.getHeaders(task) == '''
                #SBATCH -J nf-the_task_name
                #SBATCH -o /work/path/.command.log
                #SBATCH --no-requeue
                #SBATCH --signal B:USR2@30
                #SBATCH -t 00:01:00
                '''
                .stripIndent().leftTrim()

        when:
        task.config = new TaskConfig()
        task.config.time = '1h'
        task.config.memory = '50 M'
        task.config.clusterOptions = '-a 1 --signal=KILL'
        then:
        executor.getHeaders(task) == '''
                #SBATCH -J nf-the_task_name
                #SBATCH -o /work/path/.command.log
                #SBATCH --no-requeue
                #SBATCH -t 01:00:00
                #SBATCH --mem 50M
                #SBATCH -a 1 --signal=KILL
                '''
                .stripIndent().leftTrim()

        when:
        task.config = new TaskConfig()
        task.config.time = '1h'
        task.config.memory = '50 M'
        task.config.clusterOptions = ['-a 1','--signal=KILL']
        then:
        executor.getHeaders(task) == '''
                #SBATCH -J nf-the_task_name
                #SBATCH -o /work/path/.command.log
                #SBATCH --no-requeue
                #SBATCH -t 01:00:00
                #SBATCH --mem 50M
                #SBATCH -a 1
                #SBATCH --signal=KILL
                '''
            .stripIndent().leftTrim()

        when:
        task.config = new TaskConfig()
        task.config.cpus = 2
        task.config.time = '2h'
        task.config.memory = '200 M'
        task.config.clusterOptions = '-b 2'
        then:
        executor.getHeaders(task) == '''
                #SBATCH -J nf-the_task_name
                #SBATCH -o /work/path/.command.log
                #SBATCH --no-requeue
                #SBATCH --signal B:USR2@30
                #SBATCH -c 2
                #SBATCH -t 02:00:00
                #SBATCH --mem 200M
                #SBATCH -b 2
                '''
                .stripIndent().leftTrim()

        when:
        task.config = new TaskConfig()
        task.config.cpus = 8
        task.config.time = '2d 3h'
        task.config.memory = '3 G'
        task.config.clusterOptions = '-x 3'
        then:
        executor.getHeaders(task) == '''
                #SBATCH -J nf-the_task_name
                #SBATCH -o /work/path/.command.log
                #SBATCH --no-requeue
                #SBATCH --signal B:USR2@30
                #SBATCH -c 8
                #SBATCH -t 51:00:00
                #SBATCH --mem 3072M
                #SBATCH -x 3
                '''
                .stripIndent().leftTrim()

        when: 'with perCpuMemAllocation'
        executor.@perCpuMemAllocation = true
        task.config = new TaskConfig()
        task.config.cpus = 8
        task.config.memory = '24 GB'
        then:
        executor.getHeaders(task) == '''
                #SBATCH -J nf-the_task_name
                #SBATCH -o /work/path/.command.log
                #SBATCH --no-requeue
                #SBATCH --signal B:USR2@30
                #SBATCH -c 8
                #SBATCH --mem-per-cpu 3072M
                '''
                .stripIndent().leftTrim()

        when: 'with job array'
        def taskArray = Mock(TaskArrayRun) {
            config >> new TaskConfig()
            name >> task.name
            workDir >> task.workDir
            getArraySize() >> 5
        }
        then:
        executor.getHeaders(taskArray) == '''
                #SBATCH --array 0-4
                #SBATCH -J nf-the_task_name
                #SBATCH -o /dev/null
                #SBATCH --no-requeue
                #SBATCH --signal B:USR2@30
                '''
                .stripIndent().leftTrim()
    }

    def testWorkDirWithBlanks() {

        setup:
        // LSF executor
        def executor = Spy(SlurmExecutor)
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
        then:
        executor.getHeaders(task) == '''
                #SBATCH -J nf-the_task_name
                #SBATCH -o "/work/some\\ data/path/.command.log"
                #SBATCH --no-requeue
                #SBATCH --signal B:USR2@30
                '''
                .stripIndent().leftTrim()

    }


    def testQstatCommand() {

        setup:
        def executor = [:] as SlurmExecutor
        def text =
                """
                5 PD
                6 PD
                13 R
                14 CA
                15 F
                4 R
                22 S
                """.stripIndent().trim()


        when:
        def result = executor.parseQueueStatus(text)
        then:
        result.size() == 7
        result['4'] == AbstractGridExecutor.QueueStatus.RUNNING
        result['5'] == AbstractGridExecutor.QueueStatus.PENDING
        result['6'] == AbstractGridExecutor.QueueStatus.PENDING
        result['13'] == AbstractGridExecutor.QueueStatus.RUNNING
        result['14'] == AbstractGridExecutor.QueueStatus.ERROR
        result['15'] == AbstractGridExecutor.QueueStatus.ERROR
        result['22'] == AbstractGridExecutor.QueueStatus.HOLD

    }

    def testQueueStatusCommand() {
        when:
        def usr = System.getProperty('user.name')
        def exec = [:] as SlurmExecutor
        then:
        usr
        exec.queueStatusCommand(null) == ['squeue','--noheader','-o','%i %t','-t','all','-u', usr]
        exec.queueStatusCommand('xxx') == ['squeue','--noheader','-o','%i %t','-t','all','-p','xxx','-u', usr]
    }

    def 'should get array index name and start' () {
        given:
        def executor = Spy(SlurmExecutor)
        expect:
        executor.getArrayIndexName() == 'SLURM_ARRAY_TASK_ID'
        executor.getArrayIndexStart() == 0
    }

    @Unroll
    def 'should get array task id' () {
        given:
        def executor = Spy(SlurmExecutor)
        expect:
        executor.getArrayTaskId(JOB_ID, TASK_INDEX) == EXPECTED

        where:
        JOB_ID      | TASK_INDEX    | EXPECTED
        'foo'       | 1             | 'foo_1'
        'bar'       | 2             | 'bar_2'
    }

    @Unroll
    def 'should set slurm account' () {
        given:
        // task
        def task = new TaskRun()
        task.workDir = Paths.get('/work/dir')
        task.processor = Mock(TaskProcessor)
        task.processor.getSession() >> Mock(Session)
        task.config = Mock(TaskConfig)
        and:
        def executor = Spy(SlurmExecutor)
        executor.getJobNameFor(_) >> 'foo'
        executor.getName() >> 'slurm'
        executor.getSession() >> Mock(Session) { getExecConfigProp('slurm', 'account',null)>>ACCOUNT }

        when:
        def result = executor.getDirectives(task, [])
        then:
        result == EXPECTED

        where:
        ACCOUNT             | EXPECTED
        null                | ['-J', 'foo', '-o', '/work/dir/.command.log', '--no-requeue', '', '--signal', 'B:USR2@30']
        'project-123'       | ['-J', 'foo', '-o', '/work/dir/.command.log', '--no-requeue', '', '--signal', 'B:USR2@30', '-A', 'project-123']
    }
}
