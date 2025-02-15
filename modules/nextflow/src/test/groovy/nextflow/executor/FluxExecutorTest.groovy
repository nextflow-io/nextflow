/*
 * Copyright 2013-2024, Seqera Labs
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
 * @author Vanessa Sochat <sochat1@llnl.gov>
 */
class FluxExecutorTest extends Specification {

    def testParseJob() {

        given:
        def exec = [:] as FluxExecutor

        expect:
        exec.parseJobId('ƒ6nmtLVpBq') == 'ƒ6nmtLVpBq'
    }

    def testParseJobID() {

        given:
        def exec = [:] as FluxExecutor

        expect:
        exec.parseJobId('f6nmtLVpBq') == 'f6nmtLVpBq'
    }

    def testKill() {

        given:
        def executor = [:] as FluxExecutor
        expect:
        executor.killTaskCommand(123) == ['flux','job', 'cancel', '123']

    }

    def testGetCommandLine() {

        given:
        def session = Mock(Session) {
            getConfig() >> [:]
        }
        and:
        def executor = new FluxExecutor(session: session)

        // mock process
        def proc = Mock(TaskProcessor)

        // task object
        def task = new TaskRun()
        task.processor = proc
        task.name = 'my task'
        task.workDir = Paths.get('/work/path')
        task.config = new TaskConfig()

        when:
        task.config = new TaskConfig()
        task.config.queue = 'delta'
        then:
        // Flux doesn't have script headers
        executor.getHeaders(task) == ''
        executor.getSubmitCommandLine(task, Paths.get('/some/path/job.sh')) == ['flux', 'submit', '--setattr=cwd=/work/path', '--job-name="nf-my_task"', '--output=/work/path/.command.log', '--queue=delta', '/bin/bash', 'job.sh']

        when:
        task.config = new TaskConfig()
        task.config.time = '1m'
        then:
        executor.getSubmitCommandLine(task, Paths.get('/some/path/job.sh')) == ['flux', 'submit', '--setattr=cwd=/work/path', '--job-name="nf-my_task"', '--output=/work/path/.command.log', '--time-limit=01', '/bin/bash', 'job.sh']

        when:
        task.config = new TaskConfig()
        task.config.time = '1h'
        task.config.clusterOptions = '--tasks-per-node=4'
        then:
        executor.getSubmitCommandLine(task, Paths.get('/some/path/job.sh')) == ['flux', 'submit', '--setattr=cwd=/work/path', '--job-name="nf-my_task"', '--output=/work/path/.command.log', '--time-limit=60', '--tasks-per-node=4', '/bin/bash', 'job.sh']

        when:
        task.config = new TaskConfig()
        task.config.time = '1h'
        task.config.clusterOptions = '--tasks-per-node=4 --cpus-per-node=4'
        then:
        executor.getSubmitCommandLine(task, Paths.get('/some/path/job.sh')) == ['flux', 'submit', '--setattr=cwd=/work/path', '--job-name="nf-my_task"', '--output=/work/path/.command.log', '--time-limit=60', '--tasks-per-node=4', '--cpus-per-node=4', '/bin/bash', 'job.sh']

    }

    def testSubmitCommandWithTerminalOutput() {
        given:
        def session = Mock(Session) {
            getConfig() >> [flux:[terminalOutput: true]]
        }
        and:
        def executor = new FluxExecutor(session: session)
        // mock process
        def proc = Mock(TaskProcessor)
        // task object
        def task = new TaskRun()
        task.processor = proc
        task.name = 'my task'
        task.workDir = Paths.get('/work/path')
        task.config = new TaskConfig()

        expect:
        executor.getSubmitCommandLine(task, Paths.get('/some/path/job.sh')) == ['flux', 'submit', '--setattr=cwd=/work/path', '--job-name="nf-my_task"', '/bin/bash', 'job.sh']

    }

    def testWorkDirWithBlanks() {

        given:
        def session = Mock(Session) {
            getConfig() >> [:]
        }
        and:
        def executor = new FluxExecutor(session: session)
        def proc = Mock(TaskProcessor)
        def task = new TaskRun()
        task.processor = proc
        task.workDir = Paths.get('/work/some data/path')
        task.name = 'my task'

        when:
        task.index = 21
        task.config = new TaskConfig()
        then:
        executor.getSubmitCommandLine(task, Paths.get('/some/path/job.sh')) == ['flux', 'submit', '--setattr=cwd="/work/some\\ data/path"', '--job-name="nf-my_task"', '--output="/work/some\\ data/path/.command.log"', '/bin/bash', 'job.sh']

    }


    def testQstatCommand() {

        setup:
        def executor = [:] as FluxExecutor
        def text =
                """
  ƒ6upwy2MY3 R
  ƒ6upcbFjvf D
  ƒ6uon2RGVV S
  ƒ6upwy2MY4 C
  ƒ6upcbFjvh I
                """.stripIndent().trim()

        when:
        def result = executor.parseQueueStatus(text)
        then:
        result.size() == 5
        result['ƒ6upwy2MY3'] == AbstractGridExecutor.QueueStatus.RUNNING
        result['ƒ6upcbFjvf'] == AbstractGridExecutor.QueueStatus.HOLD
        result['ƒ6uon2RGVV'] == AbstractGridExecutor.QueueStatus.PENDING
        result['ƒ6upwy2MY4'] == AbstractGridExecutor.QueueStatus.DONE
        result['ƒ6upcbFjvh'] == AbstractGridExecutor.QueueStatus.DONE
    }

    def testQueueStatusCommand() {
        when:
        def usr = System.getProperty('user.name')
        def executor = [:] as FluxExecutor
        then:
        usr
        executor.queueStatusCommand(null) == ['sh', '-c', "flux jobs --suppress-header --format=\"{id.f58} {status_abbrev}\" --since=\"-15m\" --user=" + usr]
        executor.queueStatusCommand('xxx') == ['sh', '-c', "flux jobs --suppress-header --format=\"{id.f58} {status_abbrev}\" --since=\"-15m\" --queue=xxx --user=" + usr]
    }
}
