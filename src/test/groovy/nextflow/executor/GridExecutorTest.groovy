/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
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
import java.nio.file.Files

import nextflow.processor.TaskConfig
import nextflow.processor.TaskRun
import nextflow.processor.TaskStatus
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

        def executor = Mock(AbstractGridExecutor)

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

        def executor = Mock(AbstractGridExecutor)

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

        def executor = Mock(AbstractGridExecutor)

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

        def executor = Mock(AbstractGridExecutor)
        executor.checkActiveStatus(_) >> { return true }

        when:
        def handler = new GridTaskHandler(task, executor)
        handler.status = TaskStatus.RUNNING
        handler.exitFile.text = ''
        handler.exitStatusReadTimeoutMillis = 1000

        then:
        // the first try return false
        !handler.checkIfCompleted()
        // wait more the timeout defined by the property 'exitStatusReadTimeoutMillis'
        sleep 1_500
        // now 'checkIfCompleted' returns true
        handler.checkIfCompleted()
        handler.status == TaskStatus.COMPLETED
        // but the 'exitStatus' not ZERO
        handler.task.exitStatus == Integer.MAX_VALUE

    }


    def testCheckIfTerminateEmptyWithLatency() {

        setup:
        def task = new TaskRun(name: 'task1')
        task.workDir = Files.createTempDirectory('testHandler')

        def executor = Mock(AbstractGridExecutor)
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
