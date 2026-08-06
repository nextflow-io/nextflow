/*
 * Copyright 2013-2026, Seqera Labs
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

import java.nio.file.Path
import java.nio.file.Paths

import nextflow.container.DockerConfig
import nextflow.exception.ProcessFailedException
import nextflow.exception.ProcessNonZeroExitStatusException
import nextflow.file.FileHelper
import nextflow.processor.TaskBean
import nextflow.processor.TaskConfig
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import spock.lang.Specification
import test.TestHelper
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class GridTaskHandlerTest extends Specification {

    def 'should check retry predicate' () {
        given:
        def handler = new GridTaskHandler()

        when:
        def predicate = handler.retryCondition("Socket timed out")
        then:
        predicate.test(new ProcessNonZeroExitStatusException('Error', 'Socket timed out', 1, null))
        and:
        predicate.test(new ProcessNonZeroExitStatusException('Error', 'error\nBatch job submission failed\nSocket timed out on send/recv operation', 1, null))
        and:
        !predicate.test(new ProcessNonZeroExitStatusException('Error', 'OK', 0, null))

    }

    def 'should capture error cause' () {
        given:
        def task = new TaskRun(name: 'foo', workDir: Paths.get('/some/work'))
        def exec = Mock(AbstractGridExecutor) {
            getConfig() >> new ExecutorConfig([:])
        }
        def handler = Spy(new GridTaskHandler(task, exec))

        when:
        handler.submit()

        then:
        handler.safeExecute( _ )  >> { throw new ProcessNonZeroExitStatusException("Submit failed", "The limit is invalid", 10, 'qsub foo') }
        and:
        exec.createBashWrapperBuilder(task) >> Mock(BashWrapperBuilder)
        exec.pipeLauncherScript() >> false
        and:
        handler.fusionEnabled() >> false
        handler.createProcessBuilder() >> GroovyMock(ProcessBuilder)
        and:
        thrown(ProcessFailedException)
        and:
        task.stdout ==  "The limit is invalid"
        task.exitStatus == 10

    }

    def 'should replace log file with /dev/null' () {
        given:
        def WORK_DIR = Path.of('/some/dir')
        and:
        def task = Mock(TaskRun) {
            getWorkDir() >> WORK_DIR
        }
        def exec = Mock(AbstractGridExecutor) {
            getConfig() >> new ExecutorConfig([:])
        }
        def handler = Spy(new GridTaskHandler(task, exec))

        when:
        def result = handler.submitDirective(task)

        then:
        1 * exec.getHeaders(task) >> "#FOO this\n#BAR that\n#OUT file=${WORK_DIR}/.command.log\n"
        and:
        result == """\
            #FOO this
            #BAR that
            #OUT file=/dev/null
            """.stripIndent()
    }

    def 'should get stdin fusion script' () {
        given:
        def WORK_DIR = FileHelper.asPath('http://foo.com/some/dir')
        def logFile = TestHelper.createInMemTempFile('log')
        and:
        def task = Mock(TaskRun) {
            getWorkDir() >> WORK_DIR
            getLogFile() >> logFile
            getContainer() >> 'ubuntu:latest'
            getProcessor() >> Mock(TaskProcessor)
            getContainerConfig() >> Mock(DockerConfig)
            toTaskBean() >> Mock(TaskBean) { getWorkDir()>>WORK_DIR; getInputFiles()>>[:] }
            getConfig() >> Mock(TaskConfig) { getContainerOptions() >> '--this=that' }
        }
        def exec = Mock(AbstractGridExecutor) {
            getConfig() >> new ExecutorConfig([:])
        }
        def handler = Spy(new GridTaskHandler(task, exec))

        when:
        def result = handler.fusionStdinWrapper()
        then:
        handler.fusionEnabled() >> true
        exec.getHeaders(task) >> '#$ directive=one\n'
        and:
        result == '''\
                #!/bin/bash
                #$ directive=one
                docker run -i -e "FUSION_WORK=/fusion/http/foo.com/some/dir" -e "FUSION_TAGS=[.command.*|.exitcode|.fusion.*](nextflow.io/metadata=true),[*](nextflow.io/temporary=true)" --this=that ubuntu:latest /usr/bin/fusion bash '/fusion/http/foo.com/some/dir/.command.run'
                '''.stripIndent(true)
    }

    def 'should create launch command' () {
        given:
        def exec = Spy(GridTaskHandler)

        expect:
        exec.launchCmd0(new ProcessBuilder().command(['qsub', '/some/file']), null) == 'qsub /some/file'
        and:
        exec.launchCmd0(new ProcessBuilder().command(['qsub']), 'docker run /some/file') ==
                '''\
                cat << 'LAUNCH_COMMAND_EOF' | qsub
                docker run /some/file
                LAUNCH_COMMAND_EOF
                '''.stripIndent()
    }

    def 'should propagate the array index variable to a containerised array child launcher' () {
        given: 'a plain child builder built from a non-array TaskRun (arrayIndexName not set - the problematic state)'
        def childBuilder = new BashWrapperBuilder(new TaskBean(name: 'foo', workDir: Paths.get('/work/dir'), inputFiles: [:]))
        and:
        def task = Mock(TaskRun) { getWorkDir() >> Paths.get('/work/dir') }
        def exec = Mock(SlurmExecutor) {
            getConfig() >> new ExecutorConfig([:])
            getName() >> 'slurm'
            getArrayIndexName() >> 'SLURM_ARRAY_TASK_ID'
            createBashWrapperBuilder(task) >> childBuilder
        }
        def handler = Spy(new GridTaskHandler(task, exec))
        handler.withArrayChild(true)
        handler.fusionEnabled() >> false

        expect: 'the child bean does not carry the array index variable on its own'
        childBuilder.arrayIndexName == null

        when:
        def builder = handler.createTaskWrapper(task)
        then: 'the handler injects the executor array index name so the container env can expose it'
        builder.arrayIndexName == 'SLURM_ARRAY_TASK_ID'
    }

    def 'should not add the array index variable for a non-array task' () {
        given:
        def childBuilder = new BashWrapperBuilder(new TaskBean(name: 'foo', workDir: Paths.get('/work/dir'), inputFiles: [:]))
        and:
        def task = Mock(TaskRun) { getWorkDir() >> Paths.get('/work/dir') }
        def exec = Mock(SlurmExecutor) {
            getConfig() >> new ExecutorConfig([:])
            getName() >> 'slurm'
            createBashWrapperBuilder(task) >> childBuilder
        }
        def handler = Spy(new GridTaskHandler(task, exec))   // isArrayChild == false
        handler.fusionEnabled() >> false

        when:
        def builder = handler.createTaskWrapper(task)
        then:
        builder.arrayIndexName == null
        and: 'the executor array index name is never requested for a non-array task'
        0 * exec.getArrayIndexName()
    }
}
