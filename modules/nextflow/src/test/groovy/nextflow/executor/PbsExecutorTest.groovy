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

import nextflow.processor.TaskArrayRun
import nextflow.Session
import nextflow.processor.TaskConfig
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class PbsExecutorTest extends Specification {

    def testGetCommandLine() {

        given:
        def executor = Spy(PbsExecutor)
        def task = Mock(TaskRun); task.getName() >> 'hello world'

        expect:
        executor.getSubmitCommandLine(task, Paths.get('/some/path/script.sh') ) == ['qsub', '-N', 'nf-hello_world', 'script.sh']

    }

    def 'test job script headers'() {

        setup:
        def executor = Spy(PbsExecutor)
        executor.getSession() >> Mock(Session)

        // mock process
        def proc = Mock(TaskProcessor)
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
                #PBS -j oe
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
                #PBS -j oe
                #PBS -q alpha
                #PBS -l walltime=00:01:00
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
                #PBS -j oe
                #PBS -q alpha
                #PBS -l walltime=00:01:00
                #PBS -l mem=1mb
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
                #PBS -j oe
                #PBS -q delta
                #PBS -l nodes=1:ppn=2
                #PBS -l walltime=00:10:00
                #PBS -l mem=5mb
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
                #PBS -j oe
                #PBS -q delta
                #PBS -l nodes=1:ppn=8
                #PBS -l walltime=24:00:00
                #PBS -l mem=1gb
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
                #PBS -j oe
                #PBS -q delta
                #PBS -l walltime=54:10:00
                #PBS -l mem=2gb
                '''
                .stripIndent().leftTrim()

        when:
        task.config = new TaskConfig()
        task.config.clusterOptions = '-l nodes=1:x86:ppn=4'
        task.config.cpus = 2
        task.config.memory = '8g'
        then:
        executor.getHeaders(task) == '''
                #PBS -N nf-task_name
                #PBS -o /work/dir/.command.log
                #PBS -j oe
                #PBS -l mem=8gb
                #PBS -l nodes=1:x86:ppn=4
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
                #PBS -t 0-4
                #PBS -N nf-task_name
                #PBS -o /dev/null
                #PBS -j oe
                '''
                .stripIndent().leftTrim()
    }

    def WorkDirWithBlanks() {

        setup:
        def executor = Spy(PbsExecutor)
        executor.getSession() >> Mock(Session)

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
                #PBS -j oe
                '''
                .stripIndent().leftTrim()

    }

    @Unroll
    def 'should return valid job name given #name'() {
        given:
        def executor = [:] as PbsExecutor
        def task = Mock(TaskRun)
        task.getName() >> name

        expect:
        executor.getJobNameFor(task) == expected
        executor.getJobNameFor(task).size() <= 15

        where:
        name        | expected
        'hello'     | 'nf-hello'
        '12 45'     | 'nf-12_45'
        'hello(123)'| 'nf-hello123'
        'very-long-task-name-taking-more-than-15-chars' | 'nf-very-long-ta'

    }

    def testParseJobId() {

        given:
        def executor = [:] as PbsExecutor

        expect:
        executor.parseJobId('\n10.localhost\n') == '10.localhost'
        executor.parseJobId('1584288.biocluster.igb.illinois.edu') == '1584288.biocluster.igb.illinois.edu'

        when:
        executor.parseJobId('foo\nbar\n1584288.biocluster.igb.illinois.edu')
        then:
        thrown(IllegalArgumentException)


    }


    def testKillTaskCommand() {

        given:
        def executor = [:] as PbsExecutor
        expect:
        executor.killTaskCommand('100.localhost') == ['qdel', '100.localhost']

    }

    def testParseQueueStatus() {

        setup:
        def executor = [:] as PbsExecutor
        def text =
                """
                Job Id: 12.localhost
                    job_state = C
                Job Id: 13.localhost
                    job_state = R
                Job Id: 14.localhost
                    job_state = Q
                Job Id: 15.localhost
                    job_state = S
                Job Id: 16.localhost
                    job_state = E
                Job Id: 17.localhost
                    job_state = H

                """.stripIndent().trim()

        when:
        def result = executor.parseQueueStatus(text)
        then:
        result.size() == 6
        result['12.localhost'] == AbstractGridExecutor.QueueStatus.DONE
        result['13.localhost'] == AbstractGridExecutor.QueueStatus.RUNNING
        result['14.localhost'] == AbstractGridExecutor.QueueStatus.PENDING
        result['15.localhost'] == AbstractGridExecutor.QueueStatus.HOLD
        result['16.localhost'] == AbstractGridExecutor.QueueStatus.UNKNOWN
        result['17.localhost'] == AbstractGridExecutor.QueueStatus.HOLD

    }

    def 'should fetch the value' () {
        expect:
        PbsExecutor.fetchValue('Job Id:', 'Job Id:1234') == '1234'
        PbsExecutor.fetchValue('Job Id:', 'Job Id: 1234 ') == '1234'
        PbsExecutor.fetchValue('Job Id:', '  Job Id:  1234') == '1234'
    }

    def 'should return qstat command line' () {
        given:
        def executor = [:] as PbsExecutor

        expect:
        executor.queueStatusCommand(null) == ['bash','-c', "set -o pipefail; qstat -f -1 | { grep -E '(Job Id:|job_state =)' || true; }"]
        executor.queueStatusCommand('xxx') == ['bash','-c', "set -o pipefail; qstat -f -1 xxx | { grep -E '(Job Id:|job_state =)' || true; }"]
        executor.queueStatusCommand('xxx').each { assert it instanceof String }
    }

    def 'should match cluster options' () {
        expect:
        PbsExecutor.matchOptions('-l foo')
        PbsExecutor.matchOptions('-lfoo')
        PbsExecutor.matchOptions('-x -l foo')
        PbsExecutor.matchOptions('-x -lfoo')
        and:
        !PbsExecutor.matchOptions(null)
        !PbsExecutor.matchOptions('')
        !PbsExecutor.matchOptions('-x-l foo')
    }

    def 'should get array index name and start' () {
        given:
        def executor = Spy(PbsExecutor)
        expect:
        executor.getArrayIndexName() == 'PBS_ARRAYID'
        executor.getArrayIndexStart() == 0
    }

    @Unroll
    def 'should get array task id' () {
        given:
        def executor = Spy(PbsExecutor)
        expect:
        executor.getArrayTaskId(JOB_ID, TASK_INDEX) == EXPECTED

        where:
        JOB_ID      | TASK_INDEX    | EXPECTED
        'foo[]'     | 1             | 'foo[1]'
        'bar[]'     | 2             | 'bar[2]'
    }
    
    def 'should set pbs account' () {
        given:
        // task
        def task = new TaskRun()
        task.workDir = Paths.get('/work/dir')
        task.processor = Mock(TaskProcessor)
        task.processor.getSession() >> Mock(Session)
        task.config = Mock(TaskConfig)  { getClusterOptionsAsList()>>[] }
        and:
        def executor = Spy(PbsExecutor)
        executor.getJobNameFor(_) >> 'foo'
        executor.getName() >> 'pbs'
        executor.getSession() >> Mock(Session) { getExecConfigProp('pbs', 'account',null)>>ACCOUNT }

        when:
        def result = executor.getDirectives(task, [])
        then:
        result == EXPECTED

        where:
        ACCOUNT             | EXPECTED
        null                | ['-N', 'foo', '-o', '/work/dir/.command.log', '-j', 'oe']
        'project-123'       | ['-N', 'foo', '-o', '/work/dir/.command.log', '-j', 'oe', '-P', 'project-123']
    }
}
