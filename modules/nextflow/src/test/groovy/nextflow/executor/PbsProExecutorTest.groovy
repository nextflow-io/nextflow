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

import spock.lang.Specification

import java.nio.file.Paths

import nextflow.processor.TaskConfig
import nextflow.processor.TaskRun
/**
 *
 * @author Lorenz Gerber <lorenzottogerber@gmail.com>
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class PbsProExecutorTest extends Specification {

    def 'should get directives' () {
        given:
        def executor = Spy(PbsProExecutor)
        def WORK_DIR = Paths.get('/here')

        def task = Mock(TaskRun)
        task.getWorkDir() >> WORK_DIR
        task.getConfig() >> new TaskConfig()

        when:
        def result = executor.getDirectives(task, [])
        then:
        1 * executor.getJobNameFor(task) >> 'my-name'
        1 * executor.quote( WORK_DIR.resolve(TaskRun.CMD_LOG)) >> '/here/.command.log'

        result == [
                '-N', 'my-name',
                '-o', '/here/.command.log',
                '-j', 'oe'
        ]
    }

    def 'should get directives with queue' () {
        given:
        def executor = Spy(PbsProExecutor)
        def WORK_DIR = Paths.get('/foo/bar')

        def task = Mock(TaskRun)
        task.getWorkDir() >> WORK_DIR
        task.getConfig() >> new TaskConfig([ queue: 'my-queue' ])

        when:
        def result = executor.getDirectives(task, [])
        then:
        1 * executor.getJobNameFor(task) >> 'foo'
        1 * executor.quote( WORK_DIR.resolve(TaskRun.CMD_LOG)) >> '/foo/bar/.command.log'

        result == [
                '-N', 'foo',
                '-o', '/foo/bar/.command.log',
                '-j', 'oe',
                '-q', 'my-queue'
        ]
    }

    def 'should get directives with cpus' () {
        given:
        def executor = Spy(PbsProExecutor)
        def WORK_DIR = Paths.get('/foo/bar')

        def task = Mock(TaskRun)
        task.getWorkDir() >> WORK_DIR
        task.getConfig() >> new TaskConfig([ queue: 'my-queue', cpus:4 ])

        when:
        def result = executor.getDirectives(task, [])
        then:
        1 * executor.getJobNameFor(task) >> 'foo'
        1 * executor.quote( WORK_DIR.resolve(TaskRun.CMD_LOG)) >> '/foo/bar/.command.log'

        result == [
                '-N', 'foo',
                '-o', '/foo/bar/.command.log',
                '-j', 'oe',
                '-q', 'my-queue',
                '-l', 'select=1:ncpus=4'
        ]
    }

    def 'should get directives with mem' () {
        given:
        def executor = Spy(PbsProExecutor)
        def WORK_DIR = Paths.get('/foo/bar')

        def task = Mock(TaskRun)
        task.getWorkDir() >> WORK_DIR
        task.getConfig() >> new TaskConfig([ queue: 'my-queue', memory:'2 GB' ])

        when:
        def result = executor.getDirectives(task, [])
        then:
        1 * executor.getJobNameFor(task) >> 'foo'
        1 * executor.quote( WORK_DIR.resolve(TaskRun.CMD_LOG)) >> '/foo/bar/.command.log'

        result == [
                '-N', 'foo',
                '-o', '/foo/bar/.command.log',
                '-j', 'oe',
                '-q', 'my-queue',
                '-l', 'select=1:ncpus=1:mem=2048mb'
        ]
    }

    def 'should get directives with cpus and mem' () {
        given:
        def executor = Spy(PbsProExecutor)
        def WORK_DIR = Paths.get('/foo/bar')

        def task = Mock(TaskRun)
        task.getWorkDir() >> WORK_DIR
        task.getConfig() >> new TaskConfig([ queue: 'my-queue', memory:'1 GB', cpus: 8 ])

        when:
        def result = executor.getDirectives(task, [])
        then:
        1 * executor.getJobNameFor(task) >> 'foo'
        1 * executor.quote( WORK_DIR.resolve(TaskRun.CMD_LOG)) >> '/foo/bar/.command.log'

        result == [
                '-N', 'foo',
                '-o', '/foo/bar/.command.log',
                '-j', 'oe',
                '-q', 'my-queue',
                '-l', 'select=1:ncpus=8:mem=1024mb'
        ]
    }

    def 'should return qstat command line' () {
        given:
        def executor = [:] as PbsProExecutor

        expect:
        executor.queueStatusCommand(null) == ['bash','-c', "set -o pipefail; qstat -f | { egrep '(Job Id:|job_state =)' || true; }"]
        executor.queueStatusCommand('xxx') == ['bash','-c', "set -o pipefail; qstat -f xxx | { egrep '(Job Id:|job_state =)' || true; }"]
        executor.queueStatusCommand('xxx').each { assert it instanceof String }
    }

    def 'should parse queue status'() {

        setup:
        def executor = [:] as PbsProExecutor
        def text =
                """
                Job Id: 12.localhost
                    job_state = F
                Job Id: 13.localhost
                    job_state = R
                Job Id: 14.localhost
                    job_state = Q
                Job Id: 15.localhost
                    job_state = S
                Job Id: 16.localhost
                    job_state = E
                Job Id: 17.localhost
                    job_state = H

                """.stripIndent().trim()

        when:
        def result = executor.parseQueueStatus(text)
        then:
        result.size() == 6
        result['12.localhost'] == AbstractGridExecutor.QueueStatus.DONE
        result['13.localhost'] == AbstractGridExecutor.QueueStatus.RUNNING
        result['14.localhost'] == AbstractGridExecutor.QueueStatus.PENDING
        result['15.localhost'] == AbstractGridExecutor.QueueStatus.HOLD
        result['16.localhost'] == AbstractGridExecutor.QueueStatus.RUNNING
        result['17.localhost'] == AbstractGridExecutor.QueueStatus.HOLD

    }


}
