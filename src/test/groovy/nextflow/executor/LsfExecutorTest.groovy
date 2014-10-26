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
import java.nio.file.Files
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
class LsfExecutorTest extends Specification {


    def 'test bsub cmd line' () {

        setup:
        def folder = Files.createTempDirectory('test')
        // mock process
        def proc = Mock(TaskProcessor)
        def base = Mock(BaseScript)
        def config = new TaskConfig(base)
        // LSF executor
        def executor = [:] as LsfExecutor
        executor.taskConfig = config

        when:
        // process name
        proc.getName() >> 'task'
        // the script
        def script = folder.resolve('job.sh'); script.text = 'some content'
        // config
        config.queue = test_queue
        config.clusterOptions = "-x 1"
        config.cpus = test_cpu
        config.time = test_time
        config.memory = test_mem
        // task object
        def task = new TaskRun()
        task.processor = proc
        task.workDir = Paths.get('/xxx')
        task.index = 1

        then:
        executor.getSubmitCommandLine(task, script) == expected
        script.canExecute()

        cleanup:
        folder?.deleteDir()

        where:
        test_cpu    | test_time    | test_mem   | test_queue || expected
        null        | null         | null       | 'alpha'    || ['bsub','-cwd','/xxx','-o','/dev/null','-q', 'alpha', '-J', 'nf-task_1', '-x', '1', './job.sh']
        1           | null         | null       | 'alpha'    || ['bsub','-cwd','/xxx','-o','/dev/null','-q', 'alpha', '-n', '1', '-R', 'span[hosts=1]', '-J', 'nf-task_1', '-x', '1', './job.sh']
        1           | '1min'       | '10 MB'    | 'alpha'    || ['bsub','-cwd','/xxx','-o','/dev/null','-q', 'alpha', '-n', '1', '-R', 'span[hosts=1]', '-W', '00:01', '-M', '10', '-J', 'nf-task_1', '-x', '1', './job.sh']
        1           | '4h'         | '200 MB'   | 'gamma'    || ['bsub','-cwd','/xxx','-o','/dev/null','-q', 'gamma', '-n', '1', '-R', 'span[hosts=1]', '-W', '04:00', '-M', '200', '-J', 'nf-task_1', '-x', '1', './job.sh']
        4           | null         | '2 GB'     | 'gamma'    || ['bsub','-cwd','/xxx','-o','/dev/null','-q', 'gamma', '-n', '4', '-R', 'span[hosts=1]', '-M', '512', '-J', 'nf-task_1', '-x', '1', './job.sh']
        4           | '1d'         | '2 GB'     | 'gamma'    || ['bsub','-cwd','/xxx','-o','/dev/null','-q', 'gamma', '-n', '4', '-R', 'span[hosts=1]', '-W', '24:00', '-M', '512', '-J', 'nf-task_1', '-x', '1', './job.sh']
        8           | '2d'         | '2 GB'     | 'delta'    || ['bsub','-cwd','/xxx','-o','/dev/null','-q', 'delta', '-n', '8', '-R', 'span[hosts=1]', '-W', '48:00', '-M', '256', '-J', 'nf-task_1', '-x', '1', './job.sh']
        null        | '2d 12h 5m'  | '2 GB'     | 'delta'    || ['bsub','-cwd','/xxx','-o','/dev/null','-q', 'delta', '-W', '60:05', '-M', '2048', '-J', 'nf-task_1', '-x', '1', './job.sh']

    }


    def testParseJobId() {

        when:
        // executor stub object
        def executor = [:] as LsfExecutor
        then:
        executor.parseJobId( 'Job <2329803> is submitted to default queue <research-rh6>.' ) == '2329803'

    }

    def testKillCommand() {
        when:
        // executor stub object
        def executor = [:] as LsfExecutor
        then:
        executor.killTaskCommand('12345').join(' ') == 'bkill 12345'

    }

    def testQstatCommand() {

        setup:
        def executor = [:] as LsfExecutor
        def text =
                """
                6795348,RUN,Feb 17 13:26
                6795349,RUN,Feb 17 13:26
                6795351,PEND,Feb 17 13:26
                6795353,PSUSP,Feb 17 13:26
                6795354,EXIT,Feb 17 13:26
                """.stripIndent().trim()


        when:
        def result = executor.parseQueueStatus(text)
        then:
        result.size() == 5
        result['6795348'] == AbstractGridExecutor.QueueStatus.RUNNING
        result['6795349'] == AbstractGridExecutor.QueueStatus.RUNNING
        result['6795351'] == AbstractGridExecutor.QueueStatus.PENDING
        result['6795353'] == AbstractGridExecutor.QueueStatus.HOLD
        result['6795354'] == AbstractGridExecutor.QueueStatus.ERROR

    }


    def testQueueStatusCommand() {

        setup:
        def executor = [:] as LsfExecutor

        expect:
        executor.queueStatusCommand(null) == ['bjobs', '-o',  'JOBID STAT SUBMIT_TIME delimiter=\',\'', '-noheader']
        executor.queueStatusCommand('long') == ['bjobs', '-o',  'JOBID STAT SUBMIT_TIME delimiter=\',\'', '-noheader', '-q', 'long']

    }

}
