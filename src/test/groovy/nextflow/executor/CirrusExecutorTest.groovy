/*
 * Copyright (c) 2013-2015, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2015, Paolo Di Tommaso and the respective authors.
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

import nextflow.file.FileHelper
import nextflow.processor.TaskConfig
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CirrusExecutorTest extends Specification {

    def 'test cirrus directives'() {

        setup:
        def result
        def config
        // processor
        def proc = Mock(TaskProcessor)
        proc.getName() >> 'task123'
        // task
        def task = new TaskRun()
        task.metaClass.getHashLog = { 'a6f6aa6' }
        task.processor = proc
        task.index = 3
        // executor
        def executor = [:] as CirrusExecutor
        executor.metaClass.getAwsCredentials = { ['xxx','yyy'] }

        when:
        config = task.config = new TaskConfig()
        result = executor.getDirectives(task, [])
        then:
        result == ['--tag',
                   'NAME=nf-task123_3',
                    '--tag',
                    'uuid=a6f6aa6',
                   '--no-env',
                   '-v',
                   'AWS_ACCESS_KEY_ID=xxx',
                   '-v',
                   'AWS_SECRET_ACCESS_KEY=yyy'
        ]

        when:
        config = task.config = new TaskConfig()
        config.cpus = 2
        config.memory = '1G'
        config.queue = 'short'
        result = executor.getDirectives(task, [])
        then:
        result == ['--tag',
                   'NAME=nf-task123_3',
                   '--tag',
                   'uuid=a6f6aa6',
                   '--no-env',
                   '-v',
                   'AWS_ACCESS_KEY_ID=xxx',
                   '-v',
                   'AWS_SECRET_ACCESS_KEY=yyy',
                   '-q',
                   'short',
                   '-c',
                   '2',
                   '-m',
                   '1024']

        when:
        config = task.config = new TaskConfig()
        config.cpus = 4
        config.memory = '5G'
        config.queue = 'long'
        config.clusterOptions = "--any-other-option 'x and y'"
        result = executor.getDirectives(task, [])
        then:
        result == ['--tag',
                   'NAME=nf-task123_3',
                   '--tag',
                   'uuid=a6f6aa6',
                   '--no-env',
                   '-v',
                   'AWS_ACCESS_KEY_ID=xxx',
                   '-v',
                   'AWS_SECRET_ACCESS_KEY=yyy',
                   '-q',
                   'long',
                   '-c',
                   '4',
                   '-m',
                   '5120',
                   '--any-other-option',
                   'x and y']

    }

    def 'test cirrus ksub command line'() {

        given:
        // processor
        def proc = Mock(TaskProcessor)
        proc.getName() >> 'task123'

        // task
        def task = new TaskRun()
        task.metaClass.getHashLog = { '28612781' }
        task.config = new TaskConfig(queue: 'default')
        task.processor = proc
        task.index = 3
        // executor
        def executor = [:] as CirrusExecutor
        executor.metaClass.getAwsCredentials = { ['alpha','beta'] }

        when:
        def script = FileHelper.asPath('s3://bucket/work/script.sh')
        then:
        executor.getSubmitCommandLine(task, script) == [
                'ksub',
                '--tag',
                'NAME=nf-task123_3',
                '--tag',
                'uuid=28612781',
                '--no-env',
                '-v',
                'AWS_ACCESS_KEY_ID=alpha',
                '-v',
                'AWS_SECRET_ACCESS_KEY=beta',
                '-q',
                'default',
                '-w',
                'es3 cat s3://bucket/work/script.sh > script.sh && bash script.sh'
        ]

    }

    def 'test cirrus queue status command'() {

        setup:
        def executor = [:] as CirrusExecutor

        expect:
        executor.queueStatusCommand(null) == ['kqueue']
        executor.queueStatusCommand('long') == ['kqueue','-q','long']

    }

    def 'test cirrus kill command'() {

        when:
        def executor = [:] as CirrusExecutor
        then:
        executor.killTaskCommand(123) == ['kancel', '123']

    }

    def testParseQueueDump() {

        setup:
        def executor = [:] as CirrusExecutor
        def text =
                """
                Task	Attempt	Node	Code	State	ExpKb	PeakKb	MaxKb	CPUs	Queue	CmdLine	Tags	Priority	TimeHint	DiskExpKb	Dependencies
                87	56	23	-1	DISPATCHED	1048576	0	2048000	1	test-mta	["/bin/bash", "-c", "es3 cat s3://cbcrg-eu/work/b2/99242c8b37d10aefbdc4197216bb0e/.command.run > .command.run && bash .command.run"]	INSTANCE_TYPE=c3.large;name=nf-align_tree_2;uuid=b2/99242c	0	0	10485760	[]
                88	56	23	-1	RUNNING	1048576	0	2048000	1	test-mta	["/bin/bash", "-c", "es3 cat s3://cbcrg-eu/work/b2/99242c8b37d10aefbdc4197216bb0e/.command.run > .command.run && bash .command.run"]	INSTANCE_TYPE=c3.large;name=nf-align_tree_2;uuid=b2/99242c	0	0	10485760	[]
                89	-1	-1	-1	WAITING	1048576	0	2048000	1	test-mta	["/bin/bash", "-c", "es3 cat s3://cbcrg-eu/work/79/63f7592fa83ca0037934c809ba29c3/.command.run > .command.run && bash .command.run"]	name=nf-align_tree_1;uuid=79/63f759	0	0	10485760	[]
                90	-1	-1	-1	PENDING	1048576	0	2048000	1	test-mta	["/bin/bash", "-c", "es3 cat s3://cbcrg-eu/work/71/9c90dff59e5b059df7411ea8a7a9b9/.command.run > .command.run && bash .command.run"]	name=nf-align_tree_3;uuid=71/9c90df	0	0	10485760	[]
                91	-1	-1	-1	DONE	1048576	0	2048000	1	test-mta	["/bin/bash", "-c", "es3 cat s3://cbcrg-eu/work/ea/f5834a2b0827201407a313c7883706/.command.run > .command.run && bash .command.run"]	name=nf-align_tree_4;uuid=ea/f5834a	0	0	10485760	[]

                """.stripIndent().trim()


        when:
        def result = executor.parseQueueStatus(text)
        then:
        result.size() == 5
        result['87'] == AbstractGridExecutor.QueueStatus.PENDING
        result['88'] == AbstractGridExecutor.QueueStatus.RUNNING
        result['89'] == AbstractGridExecutor.QueueStatus.PENDING
        result['90'] == AbstractGridExecutor.QueueStatus.PENDING
        result['91'] == AbstractGridExecutor.QueueStatus.DONE


    }
}
