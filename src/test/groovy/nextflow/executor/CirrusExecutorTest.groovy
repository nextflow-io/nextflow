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
        task.processor = proc
        task.index = 3
        // executor
        def executor = [:] as CirrusExecutor

        when:
        config = task.config = new TaskConfig()
        result = executor.getDirectives(task, [])
        then:
        result == ['--tag', 'name=nf-task123_3','--no-env']

        when:
        config = task.config = new TaskConfig()
        config.cpus = 2
        config.memory = '1G'
        config.queue = 'short'
        result = executor.getDirectives(task, [])
        then:
        result == ['--tag', 'name=nf-task123_3','--no-env', '-q','short', '-c', '2', '-m', '1024']

        when:
        config = task.config = new TaskConfig()
        config.cpus = 4
        config.memory = '5G'
        config.queue = 'long'
        config.clusterOptions = "--any-other-option 'x and y'"
        result = executor.getDirectives(task, [])
        then:
        result == ['--tag', 'name=nf-task123_3','--no-env', '-q','long', '-c', '4', '-m', '5120', '--any-other-option', 'x and y']

    }

    def 'test cirrus ksub command line'() {

        given:
        // processor
        def proc = Mock(TaskProcessor)
        proc.getName() >> 'task123'
        // task
        def task = new TaskRun()
        task.config = new TaskConfig(queue: 'default')
        task.processor = proc
        task.index = 3
        // executor
        def executor = [:] as CirrusExecutor

        when:
        def script = FileHelper.asPath('s3://bucket/work/script.sh')
        then:
        executor.getSubmitCommandLine(task, script) == [
                'ksub',
                '--tag',
                'name=nf-task123_3',
                '--no-env',
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
        executor.queueStatusCommand('long') == ['kqueue']

    }

    def 'test cirrus kill command'() {

        when:
        def executor = [:] as CirrusExecutor
        then:
        executor.killTaskCommand(123) == ['kancel', '123']

    }
}
