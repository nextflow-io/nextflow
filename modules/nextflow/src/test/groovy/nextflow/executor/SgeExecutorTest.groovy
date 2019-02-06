/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

import nextflow.processor.TaskConfig
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class SgeExecutorTest extends Specification {


    def 'test qsub cmd line' () {

        given:
        def executor = [:] as SgeExecutor

        expect:
        executor.getSubmitCommandLine( Mock(TaskRun), Paths.get('/some/file/name.sh')) == ['qsub','-terse', 'name.sh']

    }

    def 'test qsub headers' () {

        given:
        def config
        // mock process
        def proc = Mock(TaskProcessor)

        def executor = [:] as SgeExecutor
        def task = new TaskRun()
        task.processor = proc
        task.workDir = Paths.get('/abc')
        task.name = 'the task name'

        when:

        // config
        config = task.config = new TaskConfig()
        config.queue = 'my-queue'
        config.name = 'task'

        then:
        executor.getHeaders(task) == '''
                #$ -wd /abc
                #$ -N nf-the_task_name
                #$ -o /abc/.command.log
                #$ -j y
                #$ -terse
                #$ -notify
                #$ -q my-queue
                '''
                .stripIndent().leftTrim()

        when:
        config = task.config = new TaskConfig()
        config.queue = 'my-queue'
        config.name = 'task'
        then:
        executor.getHeaders(task) == '''
                #$ -wd /abc
                #$ -N nf-the_task_name
                #$ -o /abc/.command.log
                #$ -j y
                #$ -terse
                #$ -notify
                #$ -q my-queue
                '''
                .stripIndent().leftTrim()


        when:
        config = task.config = new TaskConfig()
        config.queue = 'my-queue'
        config.name = 'task'
        config.time = '10s '
        config.clusterOptions = '-hard -alpha -beta'
        then:
        executor.getHeaders(task) == '''
                #$ -wd /abc
                #$ -N nf-the_task_name
                #$ -o /abc/.command.log
                #$ -j y
                #$ -terse
                #$ -notify
                #$ -q my-queue
                #$ -l h_rt=00:00:10
                #$ -hard -alpha -beta
                '''
                .stripIndent().leftTrim()



        when:
        config = task.config = new TaskConfig()
        config.queue = 'my-queue'
        config.name = 'task'
        config.time = '10m'
        config.memory = '1M'
        config.remove('clusterOptions')
        then:
        executor.getHeaders(task) == '''
                #$ -wd /abc
                #$ -N nf-the_task_name
                #$ -o /abc/.command.log
                #$ -j y
                #$ -terse
                #$ -notify
                #$ -q my-queue
                #$ -l h_rt=00:10:00
                #$ -l h_rss=1M,mem_free=1M
                '''
                .stripIndent().leftTrim()



        when:
        config = task.config = new TaskConfig()
        config.queue = 'my-queue'
        config.name = 'task'
        config.cpus = 1
        config.penv = 'smp'
        config.time = '2 m'
        config.memory = '2 M'
        then:
        executor.getHeaders(task) == '''
                #$ -wd /abc
                #$ -N nf-the_task_name
                #$ -o /abc/.command.log
                #$ -j y
                #$ -terse
                #$ -notify
                #$ -q my-queue
                #$ -pe smp 1
                #$ -l h_rt=00:02:00
                #$ -l h_rss=2M,mem_free=2M
                '''
                .stripIndent().leftTrim()

        when:
        config = task.config = new TaskConfig()
        config.queue = 'my-queue'
        config.name = 'task'
        config.cpus = 2
        config.penv = 'mpi'
        config.time = '3 d'
        config.memory = '3 g'
        then:
        executor.getHeaders(task) == '''
                #$ -wd /abc
                #$ -N nf-the_task_name
                #$ -o /abc/.command.log
                #$ -j y
                #$ -terse
                #$ -notify
                #$ -q my-queue
                #$ -pe mpi 2
                #$ -l h_rt=72:00:00
                #$ -l h_rss=3072M,mem_free=3072M
                '''
                .stripIndent().leftTrim()

        when:
        config = task.config = new TaskConfig()
        config.queue = 'my-queue'
        config.name = 'task'
        config.cpus = 4
        config.penv = 'orte'
        config.time = '1d3h'
        config.memory = '4 GB '
        then:
        executor.getHeaders(task) == '''
                #$ -wd /abc
                #$ -N nf-the_task_name
                #$ -o /abc/.command.log
                #$ -j y
                #$ -terse
                #$ -notify
                #$ -q my-queue
                #$ -pe orte 4
                #$ -l h_rt=27:00:00
                #$ -l h_rss=4096M,mem_free=4096M
                '''
                .stripIndent().leftTrim()

    }

    def testWorkDirWithBlanks() {

        given:
        def config
        // mock process
        def proc = Mock(TaskProcessor)

        def executor = [:] as SgeExecutor
        def task = new TaskRun()
        task.processor = proc
        task.workDir = Paths.get('/work/dir with/blanks')
        task.name = 'the task name'

        when:

        // config
        config = task.config = new TaskConfig()
        config.queue = 'my-queue'
        config.name = 'task'

        then:
        executor.getHeaders(task) == '''
                #$ -wd "/work/dir with/blanks"
                #$ -N nf-the_task_name
                #$ -o "/work/dir with/blanks/.command.log"
                #$ -j y
                #$ -terse
                #$ -notify
                #$ -q my-queue
                '''
                .stripIndent().leftTrim()

    }

    def testParseJobId() {

        when:
        def executor = [:] as SgeExecutor
        def textToParse = '''
            blah blah ..
            .. blah blah
            6472
            '''
        then:
        executor.parseJobId(textToParse) == '6472'
        executor.parseJobId('Your job 1258076 ("nf-formatBlast_2") has been submitted') == '1258076'
        executor.parseJobId('\nYour job 1258076 ("nf-formatBlast_2") has been submitted\n') == '1258076'
        executor.parseJobId('Your job 341683 ("STDIN") has been submitted') == '341683'
    }

    def testKillTaskCommand() {

        when:
        def executor = [:] as SgeExecutor
        then:
        executor.killTaskCommand(123) == ['qdel', '123']

    }


    def testParseQueueStatus() {

        setup:
        def executor = [:] as SgeExecutor
        def text =
        """
        job-ID  prior   name       user         state submit/start at     queue                          slots ja-task-ID
        -----------------------------------------------------------------------------------------------------------------
        7548318 0.00050 nf-exonera pditommaso   r     02/10/2014 12:30:51 long@node-hp0214.linux.crg.es      1
        7548348 0.00050 nf-exonera pditommaso   r     02/10/2014 12:32:43 long@node-hp0204.linux.crg.es      1
        7548349 0.00050 nf-exonera pditommaso   hqw   02/10/2014 12:32:56 long@node-hp0303.linux.crg.es      1
        7548904 0.00050 nf-exonera pditommaso   qw    02/10/2014 13:07:09                                    1
        7548960 0.00050 nf-exonera pditommaso   Eqw   02/10/2014 13:08:11                                    1
        """.stripIndent().trim()


        when:
        def result = executor.parseQueueStatus(text)
        then:
        result.size() == 5
        result['7548318'] == AbstractGridExecutor.QueueStatus.RUNNING
        result['7548348'] == AbstractGridExecutor.QueueStatus.RUNNING
        result['7548349'] == AbstractGridExecutor.QueueStatus.HOLD
        result['7548904'] == AbstractGridExecutor.QueueStatus.PENDING
        result['7548960'] == AbstractGridExecutor.QueueStatus.ERROR

    }

    def testParseQueueStatusFromUniva() {
        setup:
        def executor = [:] as SgeExecutor
        def text =
        """
        job-ID     prior   name       user         state submit/start at     queue                          jclass                         slots ja-task-ID
        ------------------------------------------------------------------------------------------------------------------------------------------------
              1220 1.00000 oliver-tes abria        r     08/29/2014 10:17:21 long@node-hp0115.linux.crg.es                                     8
              1258 0.17254 mouse.4689 epalumbo     r     08/29/2014 11:13:55 long@node-hp0515.linux.crg.es                                    16
              1261 0.17254 run_mappin epalumbo     qw    08/29/2014 11:28:11 short@node-ib0208bi.linux.crg.                                    4
              1262 0.17254 run_mappin epalumbo     Eqw   08/29/2014 11:28:31 short@node-ib0209bi.linux.crg.                                    4
        """.stripIndent().trim()


        when:
        def result = executor.parseQueueStatus(text)
        then:
        result.size() == 4
        result['1220'] == AbstractGridExecutor.QueueStatus.RUNNING
        result['1258'] == AbstractGridExecutor.QueueStatus.RUNNING
        result['1261'] == AbstractGridExecutor.QueueStatus.PENDING
        result['1262'] == AbstractGridExecutor.QueueStatus.ERROR
    }

    def testParseQueueDump() {

        setup:
        def executor = [:] as SgeExecutor
        def text =
                """
        job-ID  prior   name       user         state submit/start at     queue                          slots ja-task-ID
        -----------------------------------------------------------------------------------------------------------------
        7548318 0.00050 nf-exonera pditommaso   r     02/10/2014 12:30:51 long@node-hp0214.linux.crg.es      1
        7548348 0.00050 nf-exonera pditommaso   r     02/10/2014 12:32:43 long@node-hp0204.linux.crg.es      1
        7548349 0.00050 nf-exonera pditommaso   hqw   02/10/2014 12:32:56 long@node-hp0303.linux.crg.es      1
        7548904 0.00050 nf-exonera pditommaso   qw    02/10/2014 13:07:09                                    1
        7548960 0.00050 nf-exonera pditommaso   Eqw   02/10/2014 13:08:11                                    1
        """.stripIndent().trim()


        when:
        def status = executor.parseQueueStatus(text)
        then:
        executor.dumpQueueStatus(status).readLines().sort() == [
                '  job: 7548318: RUNNING',
                '  job: 7548348: RUNNING',
                '  job: 7548349: HOLD',
                '  job: 7548904: PENDING',
                '  job: 7548960: ERROR'
        ]


    }

    def testQueueStatusCommand() {

        setup:
        def executor = [:] as SgeExecutor

        expect:
        executor.queueStatusCommand(null) == ['qstat']
        executor.queueStatusCommand('long') == ['qstat','-q','long']

    }

}
