/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
 */

package nextflow.executor

import nextflow.processor.TaskConfig
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import spock.lang.Specification
import spock.lang.Unroll

import java.nio.file.Paths

/**
 *
 * @author Lorenz Gerber <lorenzottogerber@gmail.com>
 */
class PbsProExecutorTest extends Specification {

    def testGetCommandLine() {

        given:
        def executor = [:] as PbsProExecutor
        def task = Mock(TaskRun); task.getName() >> 'hello world'

        expect:
        executor.getSubmitCommandLine(task, Paths.get('/some/path/script.sh') ) == ['qsub', '-N', 'nf-hello_world', 'script.sh']

    }

    def testHeaders() {

        setup:
        def executor = [:] as PbsProExecutor

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
                #PBS -j oe
                cd /work/dir
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
                cd /work/dir
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
                #PBS -l select=mem=1mb
                #PBS -l walltime=00:01:00
                cd /work/dir
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
                #PBS -l select=1:ncpus=2:mem=5mb
                #PBS -l walltime=00:10:00
                cd /work/dir
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
                #PBS -l select=1:ncpus=8:mem=1gb
                #PBS -l walltime=24:00:00
                cd /work/dir
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
                #PBS -l select=mem=2gb
                #PBS -l walltime=54:10:00
                cd /work/dir
                '''
                .stripIndent().leftTrim()

    }

    def WorkDirWithBlanks() {

        setup:
        def executor = Spy(PbsProExecutor)

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
                cd "/work/dir\\ 1"
                '''
                .stripIndent().leftTrim()

    }

    @Unroll
    def 'should return valid job name given #name'() {
        given:
        def executor = [:] as PbsProExecutor
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
        def executor = [:] as PbsProExecutor

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
        def executor = [:] as PbsProExecutor
        expect:
        executor.killTaskCommand('100.localhost') == ['qdel', '100.localhost']

    }

    def testParseQueueStatus() {

        setup:
        def executor = [:] as PbsProExecutor
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
        PbsProExecutor.fetchValue('Job Id:', 'Job Id:1234') == '1234'
        PbsProExecutor.fetchValue('Job Id:', 'Job Id: 1234 ') == '1234'
        PbsProExecutor.fetchValue('Job Id:', '  Job Id:  1234') == '1234'
    }

    def 'should return qstat command line' () {
        given:
        def executor = [:] as PbsProExecutor

        expect:
        executor.queueStatusCommand(null) == ['sh','-c', "qstat -f | egrep '(Job Id:|job_state =)'"]
        executor.queueStatusCommand('xxx') == ['sh','-c', "qstat -f xxx | egrep '(Job Id:|job_state =)'"]
        executor.queueStatusCommand('xxx').each { assert it instanceof String }
    }

}
