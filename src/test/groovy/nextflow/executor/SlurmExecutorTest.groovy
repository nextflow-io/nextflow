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

    def testGetSubmitCmdLine() {

        given:
        // mock process
        def proc = Mock(TaskProcessor)
        proc.getName() >> { 'myJob' }
        // config object
        def script = Paths.get('/some/script.sh')
        // task object
        def task = new TaskRun(workDir: Paths.get('/work/path'), index: 33, processor: proc)
        // SLURM executor
        def exec = [:] as SlurmExecutor
        def config = task.config = new TaskConfig()

        when:
        config.cpus = test_cpus
        config.time = test_time
        config.memory = test_mem
        config.clusterOptions = test_opts
        then:
        exec.getSubmitCommandLine(task,script).join(' ') == expected

        where:
        test_cpus   | test_mem  | test_time | test_opts || expected
        null        | null      | null      | null      || 'sbatch -D /work/path -J nf-myJob_33 -o /dev/null script.sh'
        null        | null      | '1m'      | null      || 'sbatch -D /work/path -J nf-myJob_33 -o /dev/null -t 00:01:00 script.sh'
        1           | '50 M'    | '1h'      | '-a 1'    || 'sbatch -D /work/path -J nf-myJob_33 -o /dev/null -t 01:00:00 --mem 50 -a 1 script.sh'
        2           | '200 M'   | '2h'      | '-b 2'    || 'sbatch -D /work/path -J nf-myJob_33 -o /dev/null -c 2 -t 02:00:00 --mem 200 -b 2 script.sh'
        8           | '3 G'     | '2d'      | '-x 3'    || 'sbatch -D /work/path -J nf-myJob_33 -o /dev/null -c 8 -t 48:00:00 --mem 3072 -x 3 script.sh'
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
