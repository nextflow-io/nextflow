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
import java.nio.file.Files

import nextflow.processor.TaskConfig
import nextflow.processor.TaskRun
import nextflow.processor.TaskStatus
import nextflow.util.Duration
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class GridExecutorTest extends Specification {


    def testCheckIfNotRunning(){

        setup:
        def work = Files.createTempDirectory('test')
        def task = new TaskRun(workDir: work, name: 'hello', config: new TaskConfig(queue: 'gamma'))

        def executor = Mock(AbstractGridExecutor) {
            getConfig() >> new ExecutorConfig([:])
        }

        when:
        def handler = new GridTaskHandler(task, executor)
        handler.status = TaskStatus.SUBMITTED
        handler.queue = 'gamma'

        then:
        !handler.checkIfRunning()
        handler.status == TaskStatus.SUBMITTED

    }

    def testCheckIfRunning(){

        setup:
        def workDir = Files.createTempDirectory('checkIfRunTest')
        workDir.resolve(TaskRun.CMD_START).text = 'yes'
        def task = Mock(TaskRun)
        task.getWorkDir() >> workDir

        def executor = Mock(AbstractGridExecutor) {
            getConfig() >> new ExecutorConfig([:])
        }

        when:
        def handler = new GridTaskHandler(task, executor)
        handler.status = TaskStatus.SUBMITTED
        then:
        handler.checkIfRunning()
        handler.status == TaskStatus.RUNNING

        cleanup:
        workDir.deleteDir()
    }



    def testCheckIfTerminated(){

        setup:
        def handler = Spy(GridTaskHandler)
        handler.status = TaskStatus.RUNNING

        when:
        def complete = handler.checkIfCompleted()
        then:
        1 * handler.isRunning() >> true
        1 * handler.readExitStatus() >> null
        1 * handler.passSanityCheck() >> true
        handler.status == TaskStatus.RUNNING
        !complete
    }


    def testCheckIfTerminatedTrue() {

        setup:
        def task = new TaskRun()
        task.workDir = Files.createTempDirectory('testHandler')

        def executor = Mock(AbstractGridExecutor) {
            getConfig() >> new ExecutorConfig([:])
        }

        when:
        def handler = new GridTaskHandler(task, executor)
        handler.status = TaskStatus.RUNNING
        handler.exitFile.text = '33'
        then:
        handler.checkIfCompleted()
        handler.status == TaskStatus.COMPLETED
        handler.task.exitStatus == 33

    }

    def testCheckIfTerminateEmptyFile() {

        given:
        def task = new TaskRun(name: 'task1')
        task.workDir = Files.createTempDirectory('testHandler')

        def executor = Mock(AbstractGridExecutor) {
            getConfig() >> new ExecutorConfig([:])
        }

        when:
        def handler = new GridTaskHandler(task, executor)
        handler.status = TaskStatus.RUNNING
        handler.exitFile.text = ''
        // inject a fast-timeout awaiter so the test doesn't need to wait 270 seconds
        handler.@exitAwaiter = new ExitStatusAwaiter(Duration.of('1sec'))

        then:
        // the first try return false (within timeout window)
        !handler.checkIfCompleted()
        // wait past the 1-second timeout
        sleep 1_500
        // now 'checkIfCompleted' returns true
        handler.checkIfCompleted()
        handler.status == TaskStatus.COMPLETED
        // exit status is Integer.MAX_VALUE when the .exitcode file was persistently empty
        handler.task.exitStatus == Integer.MAX_VALUE

    }


    def testCheckIfTerminateEmptyWithLatency() {

        setup:
        def task = new TaskRun(name: 'task1')
        task.workDir = Files.createTempDirectory('testHandler')

        def executor = Mock(AbstractGridExecutor) {
            getConfig() >> new ExecutorConfig([:])
        }
        executor.checkActiveStatus(_) >> { true }

        when:
        def handler = new GridTaskHandler(task, executor)
        handler.status = TaskStatus.RUNNING
        handler.exitFile.text = ''

        assert handler.checkIfCompleted() == false
        sleep 500
        handler.exitFile.text = '123'

        then:
        handler.checkIfCompleted()
        handler.status == TaskStatus.COMPLETED
        handler.task.exitStatus == 123

    }


}
