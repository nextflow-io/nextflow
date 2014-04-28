/*
 * Copyright (c) 2013-2014, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2014, Paolo Di Tommaso and the respective authors.
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
import nextflow.processor.TaskHandler
import nextflow.processor.TaskRun
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
        def task = Mock(TaskRun)
        task.getWorkDir() >> work

        def config = Mock(TaskConfig)
        def executor = Mock(AbstractGridExecutor)

        when:
        def handler = new GridTaskHandler(task, config, executor)
        handler.status = TaskHandler.Status.SUBMITTED
        then:
        !handler.checkIfRunning()
        handler.status == TaskHandler.Status.SUBMITTED

    }

    def testCheckIfRunning(){

        setup:
        def workDir = Files.createTempDirectory('checkIfRunTest')
        workDir.resolve(TaskRun.CMD_START).text = 'yes'
        def task = Mock(TaskRun)
        task.getWorkDir() >> workDir

        def config = Mock(TaskConfig)
        def executor = Mock(AbstractGridExecutor)

        when:
        def handler = new GridTaskHandler(task, config, executor)
        handler.status = TaskHandler.Status.SUBMITTED
        then:
        handler.checkIfRunning()
        handler.status == TaskHandler.Status.RUNNING

        cleanup:
        workDir.deleteDir()
    }



    def testCheckIfTerminated(){

        setup:
        def work = Files.createTempDirectory('test')
        def task = Mock(TaskRun)
        task.getWorkDir() >> work
        def config = Mock(TaskConfig)
        def executor = [:] as AbstractGridExecutor
        executor.queueInterval = Duration.of('1min')

        when:
        def handler = new GridTaskHandler(task, config, executor)
        handler.status = TaskHandler.Status.RUNNING
        then:
        !handler.checkIfCompleted()
        handler.status == TaskHandler.Status.RUNNING

    }


    def testCheckIfTerminatedTrue() {

        setup:
        def task = new TaskRun()
        task.workDir = Files.createTempDirectory('testHandler')

        def config = Mock(TaskConfig)
        def executor = Mock(AbstractGridExecutor)

        when:
        def handler = new GridTaskHandler(task, config, executor)
        handler.status = TaskHandler.Status.RUNNING
        handler.exitFile.text = '33'
        then:
        handler.checkIfCompleted()
        handler.status == TaskHandler.Status.COMPLETED
        handler.task.exitStatus == 33

    }

    def testCheckIfTerminateEmptyFile() {

        given:
        def task = new TaskRun()
        task.workDir = Files.createTempDirectory('testHandler')

        def config = Mock(TaskConfig)
        def executor = Mock(AbstractGridExecutor)
        executor.checkActiveStatus(_) >> { return true }

        when:
        def handler = new GridTaskHandler(task, config, executor)
        handler.status = TaskHandler.Status.RUNNING
        handler.exitFile.text = ''
        handler.exitStatusReadTimeoutMillis = 1000

        then:
        // the first try return false
        !handler.checkIfCompleted()
        // wait more the timeout defined by the property 'exitStatusReadTimeoutMillis'
        sleep 1_500
        // now 'checkIfCompleted' returns true
        handler.checkIfCompleted()
        handler.status == TaskHandler.Status.COMPLETED
        // but the 'exitStatus' not ZERO
        handler.task.exitStatus == Integer.MAX_VALUE

    }


    def testCheckIfTerminateEmptyWithLatency() {

        setup:
        def task = new TaskRun()
        task.workDir = Files.createTempDirectory('testHandler')
        def config = Mock(TaskConfig)
        def executor = Mock(AbstractGridExecutor)
        executor.checkActiveStatus(_) >> { true }

        when:
        def handler = new GridTaskHandler(task, config, executor)
        handler.status = TaskHandler.Status.RUNNING
        handler.exitFile.text = ''

        assert handler.checkIfCompleted() == false
        sleep 500
        handler.exitFile.text = '123'

        then:
        handler.checkIfCompleted()
        handler.status == TaskHandler.Status.COMPLETED
        handler.task.exitStatus == 123

    }



}
