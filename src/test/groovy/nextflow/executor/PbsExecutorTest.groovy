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
import nextflow.script.BaseScript
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class PbsExecutorTest extends Specification {

    def testQsubCommandLine () {

        given:
        // mock process
        def proc = Mock(TaskProcessor)
        def base = Mock(BaseScript)
        def config = new TaskConfig(base)
        // LSF executor
        def executor = [:] as PbsExecutor
        executor.taskConfig = config

        when:
        // process name
        proc.getName() >> 'task_x'
        // the script
        def script = Paths.get('job.sh')
        // config
        config.queue = 'my-queue'
        config.memory = test_mem
        config.time = test_time
        config.cpus = test_cpus
        config.nodes = test_nodes
        config.clusterOptions = '-extra opt'

        def task = new TaskRun()
        task.processor = proc
        task.workDir = Paths.get('/work/dir')
        task.index = 2

        then:
        executor.getSubmitCommandLine(task,script) == expected.split(' ') as List

        where:
        test_cpus| test_nodes | test_mem | test_time|| expected
        null     | null       | null     | null     || 'qsub -d /work/dir -N nf-task_x_2 -o /dev/null -e /dev/null -V -q my-queue -extra opt job.sh'
        null     | null       | null     | '1m'     || 'qsub -d /work/dir -N nf-task_x_2 -o /dev/null -e /dev/null -V -q my-queue -l walltime=00:01:00 -extra opt job.sh'
        null     | null       | '1m'     | '1m'     || 'qsub -d /work/dir -N nf-task_x_2 -o /dev/null -e /dev/null -V -q my-queue -l walltime=00:01:00 -l mem=1mb -extra opt job.sh'
        2        | null       | '5m'     | '10m'    || 'qsub -d /work/dir -N nf-task_x_2 -o /dev/null -e /dev/null -V -q my-queue -l nodes=1:ppn=2 -l walltime=00:10:00 -l mem=5mb -extra opt job.sh'
        8        | 2          | '1g'     | '1d'     || 'qsub -d /work/dir -N nf-task_x_2 -o /dev/null -e /dev/null -V -q my-queue -l nodes=2:ppn=8 -l walltime=24:00:00 -l mem=1gb -extra opt job.sh'
        null     | 4          | '2g'     | '1d'     || 'qsub -d /work/dir -N nf-task_x_2 -o /dev/null -e /dev/null -V -q my-queue -l nodes=4:ppn=1 -l walltime=24:00:00 -l mem=2gb -extra opt job.sh'


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
        result['16.localhost'] == AbstractGridExecutor.QueueStatus.UNKNWON
        result['17.localhost'] == AbstractGridExecutor.QueueStatus.HOLD

    }

}
