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
class PbsExecutorTest extends Specification {

    def testGetCommandLine() {

        when:
        def executor = [:] as PbsExecutor
        then:
        executor.getSubmitCommandLine(Mock(TaskRun), Paths.get('(/some/path/script.sh') ) == ['qsub', 'script.sh']

    }

    def testHeaders() {

        setup:
        def executor = [:] as PbsExecutor

        // mock process
        def proc = Mock(TaskProcessor)
        // process name
        proc.getName() >> 'task'

        // task object
        def task = new TaskRun()
        task.processor = proc
        task.workDir = Paths.get('/work/dir')
        task.index = 33

        when:
        task.config = new TaskConfig()
        then:
        executor.getHeaders(task) == '''
                #PBS -d /work/dir
                #PBS -N nf-task_33
                #PBS -o /dev/null
                #PBS -e /dev/null
                #PBS -V
                '''
                .stripIndent().leftTrim()


        when:
        task.config = new TaskConfig()
        task.config.queue = 'alpha'
        task.config.time = '1m'
        then:
        executor.getHeaders(task) == '''
                #PBS -d /work/dir
                #PBS -N nf-task_33
                #PBS -o /dev/null
                #PBS -e /dev/null
                #PBS -V
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
                #PBS -d /work/dir
                #PBS -N nf-task_33
                #PBS -o /dev/null
                #PBS -e /dev/null
                #PBS -V
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
                #PBS -d /work/dir
                #PBS -N nf-task_33
                #PBS -o /dev/null
                #PBS -e /dev/null
                #PBS -V
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
                #PBS -d /work/dir
                #PBS -N nf-task_33
                #PBS -o /dev/null
                #PBS -e /dev/null
                #PBS -V
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
                #PBS -d /work/dir
                #PBS -N nf-task_33
                #PBS -o /dev/null
                #PBS -e /dev/null
                #PBS -V
                #PBS -q delta
                #PBS -l walltime=54:10:00
                #PBS -l mem=2gb
                '''
                .stripIndent().leftTrim()

    }


    def testParseJobId() {

        when:
        def executor = [:] as PbsExecutor
        def textToParse = '\n10.host\n'
        then:
        executor.parseJobId(textToParse) == '10.host'
    }


    def testKillTaskCommand() {

        when:
        def executor = [:] as PbsExecutor
        then:
        executor.killTaskCommand('100.hostname') == ['qdel', '100.hostname'] as String[]

    }

    def testParseQueueStatus() {

        setup:
        def executor = [:] as PbsExecutor
        def text =
                """
        Job id                    Name             User            Time Use S Queue
        ------------------------- ---------------- --------------- -------- - -----
        12.localhost              test.sh          vagrant                0 C batch
        13.localhost              test.sh          vagrant                0 R batch
        14.localhost              test.sh          vagrant                0 Q batch
        15.localhost              test.sh          vagrant                0 S batch
        16.localhost              test.sh          vagrant                0 E batch
        17.localhost              test.sh          vagrant                0 H batch
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

}
