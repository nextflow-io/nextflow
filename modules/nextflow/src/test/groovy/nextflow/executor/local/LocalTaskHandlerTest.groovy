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

package nextflow.executor.local

import java.nio.file.Path
import java.util.concurrent.TimeUnit

import nextflow.Global
import nextflow.container.DockerConfig
import nextflow.exception.ProcessException
import nextflow.file.http.XPath
import nextflow.processor.TaskBean
import nextflow.processor.TaskConfig
import nextflow.processor.TaskRun
import nextflow.processor.TaskStatus
import nextflow.util.Duration
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class LocalTaskHandlerTest extends Specification {

    def 'should create local process builder' () {
        given:
        def task = Mock(TaskRun) {
            getWorkDir() >> Path.of('/some/work/dir')
            getConfig() >> Mock(TaskConfig)
        }
        and:
        def handler = Spy(new LocalTaskHandler(task, Mock(LocalExecutor)))

        when:
        def builder = handler.createLaunchProcessBuilder()
        then:
        handler.fusionEnabled() >> false
        and:
        builder.command() == ['/bin/bash','-ue','.command.run']
        builder.directory() == new File('/some/work/dir')
        builder.redirectErrorStream()
        builder.redirectOutput().file() == new File('/some/work/dir/.command.log')
    }

    def 'should create fusion process builder' () {
        given:
        Global.config = [:]
        def WORK_DIR = XPath.get('http://some/work/dir')
        and:
        def bean = new TaskBean(workDir: WORK_DIR, inputFiles: [:])
        and:
        def task = Mock(TaskRun) {
            getContainer() >> 'ubuntu:latest'
            getWorkDir() >> WORK_DIR
            getConfig() >> Mock(TaskConfig)
            getContainerConfig() >> new DockerConfig(enabled:true)
            toTaskBean() >> bean
        }
        def executor = Mock(LocalExecutor)
        and:
        def handler = Spy(new LocalTaskHandler(task, executor))

        when:
        def builder = handler.createLaunchProcessBuilder()
        then:
        handler.fusionEnabled() >> true
        and:
        builder.command() == ['sh','-c','docker run -i -e "FUSION_WORK=/fusion/http/some/work/dir" -e "FUSION_TAGS=[.command.*|.exitcode|.fusion.*](nextflow.io/metadata=true),[*](nextflow.io/temporary=true)" --rm --privileged ubuntu:latest /usr/bin/fusion bash \'/fusion/http/some/work/dir/.command.run\'']
        builder.directory() == null
        builder.redirectErrorStream()
        builder.redirectOutput().file()

        cleanup:
        builder?.redirectOutput()?.file()?.delete()
    }

    def 'should kill task when task exceeds time limit' () {
        given:
        def workDir = Path.of('/tmp/test-work-dir')
        def task = Mock(TaskRun) {
            getWorkDir() >> workDir
            getConfig() >> Mock(TaskConfig) {
                getTime() >> Duration.of(100)
            }
        }
        and:
        def handler = Spy(new LocalTaskHandler(task, Mock(LocalExecutor))) {
            buildTaskWrapper() >> {}
            elapsedTimeMillis() >> 200
        }
        handler.@process = Mock(Process) {
            waitFor(_ as long, _ as TimeUnit) >> true   // process has already terminated
            exitValue() >> 143  // Typical exit code for SIGTERM
        }
        handler.status = TaskStatus.RUNNING

        when:
        def completed = handler.checkIfCompleted()

        then:
        completed == true
        1 * task.setExitStatus(143)
        1 * task.setError(_ as ProcessException)
        handler.status == TaskStatus.COMPLETED
    }

    def 'should not throw when the process has not exited yet after being destroyed' () {
        given: 'a task over its time limit whose process has not been reaped yet'
        def task = Mock(TaskRun) {
            getWorkDir() >> Path.of('/tmp/test-work-dir')
            getConfig() >> Mock(TaskConfig) { getTime() >> Duration.of(100) }
        }
        and:
        def handler = Spy(new LocalTaskHandler(task, Mock(LocalExecutor))) {
            buildTaskWrapper() >> {}
            elapsedTimeMillis() >> 200
        }
        and: 'exitValue() throws, as ProcessImpl does before the child is reaped'
        handler.@process = Mock(Process) {
            waitFor(_ as long, _ as TimeUnit) >> false
            exitValue() >> { throw new IllegalThreadStateException('process hasn\'t exited') }
        }
        handler.status = TaskStatus.RUNNING

        when:
        def completed = handler.checkIfCompleted()

        then: 'the task completes instead of propagating the exception to the poll loop'
        noExceptionThrown()
        completed == true
        1 * task.setExitStatus(Integer.MAX_VALUE)
        1 * task.setError(_ as ProcessException)
        handler.status == TaskStatus.COMPLETED
    }

    def 'should read the exit status of a real process killed on timeout' () {
        given: 'a genuinely running process, destroyed the way checkIfCompleted does'
        def task = Mock(TaskRun) {
            getWorkDir() >> Path.of('/tmp/test-work-dir')
            getConfig() >> Mock(TaskConfig) { getTime() >> Duration.of(100) }
        }
        and:
        def handler = Spy(new LocalTaskHandler(task, Mock(LocalExecutor))) {
            buildTaskWrapper() >> {}
            elapsedTimeMillis() >> 200
        }
        and:
        def proc = new ProcessBuilder('bash', '-c', 'sleep 300')
                .redirectErrorStream(true)
                .redirectOutput(new File('/dev/null'))
                .start()
        handler.@process = proc
        handler.status = TaskStatus.RUNNING

        when:
        def completed = handler.checkIfCompleted()

        then: 'SIGTERM is observed rather than IllegalThreadStateException'
        noExceptionThrown()
        completed == true
        1 * task.setExitStatus(143)
        handler.status == TaskStatus.COMPLETED

        cleanup:
        proc?.destroyForcibly()
    }
}
