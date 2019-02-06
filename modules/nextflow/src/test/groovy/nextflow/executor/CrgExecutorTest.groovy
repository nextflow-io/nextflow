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

import nextflow.Session
import nextflow.processor.TaskConfig
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CrgExecutorTest extends Specification {

    def testQsubCliCommand () {

        given:
        def executor = [:] as CrgExecutor

        expect:
        executor.getSubmitCommandLine( Mock(TaskRun), Paths.get('/some/file/name.sh')) == ['qsub', '-terse','name.sh']

    }

    def testGetDirectives() {

        setup:
        def result
        def config
        // task
        def task = new TaskRun()
        task.workDir = Paths.get('/work/dir')
        task.metaClass.getHashLog = { 'a6f6aa6' }
        task.processor = Mock(TaskProcessor)
        task.processor.getSession() >> new Session()
        task.name = 'task this and that'
        // executor
        def executor = [:] as CrgExecutor

        when:
        config = task.config = new TaskConfig()
        result = executor.getDirectives(task, [])
        then:
        result == [
                '-wd',
                '/work/dir',
                '-N',
                'nf-task_this_and_that',
                '-o',
                '/work/dir/.command.log',
                '-j',
                'y',
                '-terse',
                '',
                '-notify',
                ''
        ]

        when:
        config = task.config = new TaskConfig()
        config.cpus = 4
        result = executor.getDirectives(task, [])
        then:
        result == [
                '-wd',
                '/work/dir',
                '-N',
                'nf-task_this_and_that',
                '-o',
                '/work/dir/.command.log',
                '-j',
                'y',
                '-terse',
                '',
                '-notify',
                '',
                '-pe',
                'smp 4'
        ]

        when:
        config = task.config = new TaskConfig()
        config.cpus = 8
        config.penv = 'orte'
        result = executor.getDirectives(task, [])
        then:
        result == [
                '-wd',
                '/work/dir',
                '-N',
                'nf-task_this_and_that',
                '-o',
                '/work/dir/.command.log',
                '-j',
                'y',
                '-terse',
                '',
                '-notify',
                '',
                '-pe',
                'orte 8'
        ]

        when:
        task.processor = Mock(TaskProcessor)
        task.processor.getSession() >> new Session(docker: [enabled: true])
        config = task.config = new TaskConfig()
        config.container = 'busybox'

        result = executor.getDirectives(task, [])
        then:
        result == [
                '-wd',
                '/work/dir',
                '-N',
                'nf-task_this_and_that',
                '-o',
                '/work/dir/.command.log',
                '-j',
                'y',
                '-terse',
                '',
                '-notify',
                '',
                '-binding',
                'env linear:1',
                '-soft',
                '-l docker_images=*;busybox;*'
        ]

    }

    def testGetHeaders () {

        given:
        // LSF executor
        def executor = [:] as CrgExecutor
        executor.session = new Session()

        def task = new TaskRun()
        task.processor = Mock(TaskProcessor)
        task.processor.getSession() >> new Session()
        task.processor.getName() >> 'task_x'
        task.workDir = Paths.get('/abc')
        task.name = 'mapping tag'

        when:
        // config
        task.config = new TaskConfig(
                        queue: 'short',
                        memory: '200 MB',
                        time: '1d',
                        disk: '2G'
                    )

        then:
        executor.getHeaders(task) == '''
                    #$ -wd /abc
                    #$ -N nf-mapping_tag
                    #$ -o /abc/.command.log
                    #$ -j y
                    #$ -terse
                    #$ -notify
                    #$ -q short
                    #$ -l h_rt=24:00:00
                    #$ -l h_vmem=200M,virtual_free=200M
                    #$ -l disk=2048M
                    '''
                    .stripIndent().leftTrim()

        when:
        executor.session.config.docker = [enabled: false]
        task.config = new TaskConfig(
                queue: 'short',
                memory: '4 GB',
                time: '1d',
                container: 'ubuntu'
        )

        then:
        executor.getHeaders(task) == '''
                    #$ -wd /abc
                    #$ -N nf-mapping_tag
                    #$ -o /abc/.command.log
                    #$ -j y
                    #$ -terse
                    #$ -notify
                    #$ -q short
                    #$ -l h_rt=24:00:00
                    #$ -l h_vmem=4096M,virtual_free=4096M
                    '''
                .stripIndent().leftTrim()


        when:
        task.processor = Mock(TaskProcessor)
        task.processor.getSession() >> new Session(docker: [enabled: true])
        task.config = new TaskConfig(
                queue: 'short',
                memory: '4 GB',
                time: '1d',
                container: 'ubuntu'
        )

        then:
        executor.getHeaders(task) == '''
                    #$ -wd /abc
                    #$ -N nf-mapping_tag
                    #$ -o /abc/.command.log
                    #$ -j y
                    #$ -terse
                    #$ -notify
                    #$ -q short
                    #$ -l h_rt=24:00:00
                    #$ -l h_vmem=4096M,virtual_free=4096M
                    #$ -binding env linear:1
                    #$ -soft -l docker_images=*;ubuntu;*
                    '''
                .stripIndent().leftTrim()

        when:
        task.processor = Mock(TaskProcessor)
        task.processor.getSession() >> new Session(docker: [enabled: true])
        task.config = new TaskConfig(
                memory: '3 g',
                time: '3 d',
                cpus: '2',
                penv: 'mpi',
                queue: 'long',
                container: 'busybox'
        )

        then:
        executor.getHeaders(task) == '''
                    #$ -wd /abc
                    #$ -N nf-mapping_tag
                    #$ -o /abc/.command.log
                    #$ -j y
                    #$ -terse
                    #$ -notify
                    #$ -q long
                    #$ -pe mpi 2
                    #$ -l h_rt=72:00:00
                    #$ -l h_vmem=3072M,virtual_free=3072M
                    #$ -binding env linear:2
                    #$ -R y
                    #$ -soft -l docker_images=*;busybox;*
                    '''
                .stripIndent().leftTrim()

        /*
         * repeat the same test to make sure it returns the same result
         * i.e. check it is idempotent
         */
        when:
        executor.session.config.docker = [enabled: true]
        task.config = new TaskConfig(
                memory: '3 g',
                time: '3 d',
                cpus: '4',
                penv: 'mpi',
                queue: 'long',
                container: 'busybox',
        )

        then:
        executor.getHeaders(task) == '''
                    #$ -wd /abc
                    #$ -N nf-mapping_tag
                    #$ -o /abc/.command.log
                    #$ -j y
                    #$ -terse
                    #$ -notify
                    #$ -q long
                    #$ -pe mpi 4
                    #$ -l h_rt=72:00:00
                    #$ -l h_vmem=3072M,virtual_free=3072M
                    #$ -binding env linear:4
                    #$ -R y
                    #$ -soft -l docker_images=*;busybox;*
                    '''
                    .stripIndent().leftTrim()


    }


    def testParseJobId() {

        when:
        def executor = [:] as CrgExecutor
        def textToParse = '''
            blah blah ..
            .. blah blah
            6472
            '''
        then:
        executor.parseJobId(textToParse) == '6472'
    }

    def testKillTaskCommand() {

        when:
        def executor = [:] as CrgExecutor
        then:
        executor.killTaskCommand(123) == ['qdel', '123']

    }


    def testParseQueueStatus() {

        setup:
        def executor = [:] as CrgExecutor
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
        def executor = [:] as CrgExecutor
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
        def executor = [:] as CrgExecutor
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
        def executor = [:] as CrgExecutor

        expect:
        executor.queueStatusCommand(null) == ['qstat']
        executor.queueStatusCommand('long') == ['qstat','-q','long']

    }

    def 'executor should inject `cpuset` option on docker run command' () {
        given:
        // task
        def task = new TaskRun()
        task.workDir = Paths.get('/some/dir')
        task.script = 'echo hello'
        task.processor = Mock(TaskProcessor)
        task.processor.getSession() >> new Session()
        task.processor.getProcessEnvironment() >> [:]
        task.processor.getConfig() >> [:]
        task.name = 'the-name'
        task.config = new TaskConfig()
        def executor = [:] as CrgExecutor

        when:
        def builder = executor.createBashWrapperBuilder(task)
        then:
        builder.headerScript == '''
            #$ -wd /some/dir
            #$ -N nf-the-name
            #$ -o /some/dir/.command.log
            #$ -j y
            #$ -terse
            #$ -notify
            '''
            .stripIndent().leftTrim()


        when:
        task.config.container = 'foo'
        task.processor = Mock(TaskProcessor)
        task.processor.getProcessEnvironment() >> [:]
        task.processor.getSession() >> new Session(docker: [enabled: true])
        task.processor.getConfig() >> [:]

        builder = executor.createBashWrapperBuilder(task)
        then:
        builder.headerScript == '''
            #$ -wd /some/dir
            #$ -N nf-the-name
            #$ -o /some/dir/.command.log
            #$ -j y
            #$ -terse
            #$ -notify
            #$ -binding env linear:1
            #$ -soft -l docker_images=*;foo;*

            cpuset=${cpuset:=''}
            [[ $SGE_BINDING ]] && cpuset="--cpuset-cpus $(echo $SGE_BINDING | sed 's/ /,/g')"
            '''
                .stripIndent().leftTrim()


        when:
        task.config.container = 'foo'
        task.processor = Mock(TaskProcessor)
        task.processor.getProcessEnvironment() >> [:]
        task.processor.getSession() >> new Session(docker: [enabled: true, legacy: true])
        task.processor.getConfig() >> [:]

        builder = executor.createBashWrapperBuilder(task)
        then:
        builder.headerScript == '''
            #$ -wd /some/dir
            #$ -N nf-the-name
            #$ -o /some/dir/.command.log
            #$ -j y
            #$ -terse
            #$ -notify
            #$ -binding env linear:1
            #$ -soft -l docker_images=*;foo;*

            cpuset=${cpuset:=''}
            [[ $SGE_BINDING ]] && cpuset="--cpuset $(echo $SGE_BINDING | sed 's/ /,/g')"
            '''
                .stripIndent().leftTrim()
    }
}
