/*
 * Copyright (c) 2013-2014, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2014, Paolo Di Tommaso and the respective authors.
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
import java.nio.file.Paths

import nextflow.processor.TaskConfig
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class SlurmExecutorTest extends Specification {

    def testParseJob() {

        when:
        def exec = [:] as SlurmExecutor
        then:
        exec.parseJobId('Submitted batch job 10') == '10'

    }

    def testKill() {

        given:
        def exec = [:] as SlurmExecutor
        expect:
        exec.killTaskCommand(123) == ['scancel','123']

    }

    def testGetCommandLine() {

        when:
        def exec = [:] as SlurmExecutor
        then:
        exec.getSubmitCommandLine(Mock(TaskRun), Paths.get('/some/path/job.sh')) == ['sbatch', 'job.sh']
    }

    def testGetHeaders() {

        setup:
        // LSF executor
        def executor = [:] as SlurmExecutor

        // mock process
        def proc = Mock(TaskProcessor)
        // process name
        proc.getName() >> 'task'

        // task object
        def task = new TaskRun()
        task.processor = proc
        task.workDir = Paths.get('/work/path')

        when:
        task.index = 21
        task.config = new TaskConfig()
        then:
        executor.getHeaders(task) == '''
                #SBATCH -D /work/path
                #SBATCH -J nf-task_21
                #SBATCH -o /dev/null
                '''
                .stripIndent().leftTrim()

        when:
        task.index = 22
        task.config = new TaskConfig()
        task.config.time = '1m'
        then:
        executor.getHeaders(task) == '''
                #SBATCH -D /work/path
                #SBATCH -J nf-task_22
                #SBATCH -o /dev/null
                #SBATCH -t 00:01:00
                '''
                .stripIndent().leftTrim()

        when:
        task.index = 23
        task.config = new TaskConfig()
        task.config.time = '1h'
        task.config.memory = '50 M'
        task.config.clusterOptions = '-a 1'

        then:
        executor.getHeaders(task) == '''
                #SBATCH -D /work/path
                #SBATCH -J nf-task_23
                #SBATCH -o /dev/null
                #SBATCH -t 01:00:00
                #SBATCH --mem 50
                #SBATCH -a 1
                '''
                .stripIndent().leftTrim()

        when:
        task.index = 24
        task.config = new TaskConfig()
        task.config.cpus = 2
        task.config.time = '2h'
        task.config.memory = '200 M'
        task.config.clusterOptions = '-b 2'

        then:
        executor.getHeaders(task) == '''
                #SBATCH -D /work/path
                #SBATCH -J nf-task_24
                #SBATCH -o /dev/null
                #SBATCH -c 2
                #SBATCH -t 02:00:00
                #SBATCH --mem 200
                #SBATCH -b 2
                '''
                .stripIndent().leftTrim()

        when:
        task.index = 55
        task.config = new TaskConfig()
        task.config.cpus = 8
        task.config.time = '2d 3h'
        task.config.memory = '3 G'
        task.config.clusterOptions = '-x 3'

        then:
        executor.getHeaders(task) == '''
                #SBATCH -D /work/path
                #SBATCH -J nf-task_55
                #SBATCH -o /dev/null
                #SBATCH -c 8
                #SBATCH -t 51:00:00
                #SBATCH --mem 3072
                #SBATCH -x 3
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
        def exec = [:] as SlurmExecutor
        then:
        exec.queueStatusCommand(null) == ['squeue','-h','-o','%i %t']
        exec.queueStatusCommand('xxx') == ['squeue','-h','-o','%i %t']


    }
}
