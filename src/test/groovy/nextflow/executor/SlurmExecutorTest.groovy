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

    def testGetSubmitCmdLine() {

        given:
        // mock process
        def proc = Mock(TaskProcessor)
        proc.getName() >> { 'myJob' }
        // mock script
        def base = Mock(BaseScript)
        // config object
        def config = new TaskConfig(base)
        def script = Paths.get('/some/script.sh')
        // task object
        def task = new TaskRun(workDir: Paths.get('/work/path'), index: 33, processor: proc)
        // SLURM executor
        def exec = [:] as SlurmExecutor
        exec.taskConfig = config

        when:
        config.cpus = test_cpus
        config.nodes = test_nodes
        config.time = test_time
        config.clusterOptions = '-x -y -z'
        then:
        exec.getSubmitCommandLine(task,script).join(' ') == expected

        where:
        test_cpus   | test_nodes| test_time || expected
        null        | null      | null      || 'sbatch -D /work/path -J nf-myJob_33 -o /dev/null -x -y -z script.sh'
        null        | null      | '1m'      || 'sbatch -D /work/path -J nf-myJob_33 -o /dev/null -t 00:01:00 -x -y -z script.sh'
        1           | 1         | '1h'      || 'sbatch -D /work/path -J nf-myJob_33 -o /dev/null -c 1 -N 1 -t 01:00:00 -x -y -z script.sh'
        2           | 1         | '2h'      || 'sbatch -D /work/path -J nf-myJob_33 -o /dev/null -c 2 -N 1 -t 02:00:00 -x -y -z script.sh'
        8           | 4         | '2d'      || 'sbatch -D /work/path -J nf-myJob_33 -o /dev/null -c 8 -N 4 -t 48:00:00 -x -y -z script.sh'
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
                """.stripIndent().trim()


        when:
        def result = executor.parseQueueStatus(text)
        then:
        result.size() == 6
        result['4'] == AbstractGridExecutor.QueueStatus.RUNNING
        result['5'] == AbstractGridExecutor.QueueStatus.PENDING
        result['6'] == AbstractGridExecutor.QueueStatus.PENDING
        result['13'] == AbstractGridExecutor.QueueStatus.RUNNING
        result['14'] == AbstractGridExecutor.QueueStatus.ERROR
        result['15'] == AbstractGridExecutor.QueueStatus.ERROR

    }

    def testQueueStatusCommand() {
        when:
        def exec = [:] as SlurmExecutor
        then:
        exec.queueStatusCommand(null) == ['squeue','-h','-o','%i %t']
        exec.queueStatusCommand('xxx') == ['squeue','-h','-o','%i %t']


    }
}
