/*
 * Copyright (c) 2012, the authors.
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
import nextflow.processor.TaskRun
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class LsfExecutorTest extends Specification {


    def 'test bsub cmd line' () {

        setup:
        def map = [:]
        map.queue = 'hpc-queue1'
        map.maxMemory = '2GB'
        map.maxDuration = '3h'
        map.clusterOptions = " -M 4000  -R 'rusage[mem=4000] select[mem>4000]' --X \"abc\" "
        map.name = 'task'
        def config = new TaskConfig(map)

        // dummy script
        def wrapper = Paths.get('.job.sh')
        wrapper.text = 'test script'

        // executor stub object
        def executor = [:] as LsfExecutor
        executor.taskConfig = config

        when:
        def task = new TaskRun(name: 'task 1')
        task.workDirectory = Paths.get('/xxx')

        then:
        executor.getSubmitCommandLine(task, wrapper) == ['bsub','-cwd','/xxx','-o','/dev/null','-q', 'hpc-queue1', '-J', 'nf-task_1', '-M', '4000' ,'-R' ,'rusage[mem=4000] select[mem>4000]', '--X', 'abc', './.job.sh']
        wrapper.canExecute()

        cleanup:
        wrapper?.delete()
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
