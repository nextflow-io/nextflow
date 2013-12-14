/*
 * Copyright (c) 2012, the authors.
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
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class GridExecutorTest extends Specification {


    def testCheckIfNotRunning(){

        setup:
        def task = Mock(TaskRun)
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
        def task = Mock(TaskRun)
        task.getCmdStartedFile() >> Files.createTempFile('checkIfRun',null)

        def config = Mock(TaskConfig)
        def executor = Mock(AbstractGridExecutor)

        when:
        def handler = new GridTaskHandler(task, config, executor)
        handler.status = TaskHandler.Status.SUBMITTED
        then:
        handler.checkIfRunning()
        handler.status == TaskHandler.Status.RUNNING

    }



    def testCheckIfTerminated(){

        setup:
        def task = Mock(TaskRun)
        def config = Mock(TaskConfig)
        def executor = Mock(AbstractGridExecutor)

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
        task.workDirectory = Files.createTempDirectory('testHandler')

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

        def task = new TaskRun()
        task.workDirectory = Files.createTempDirectory('testHandler')

        def config = Mock(TaskConfig)
        def executor = Mock(AbstractGridExecutor)

        when:
        def handler = new GridTaskHandler(task, config, executor)
        handler.status = TaskHandler.Status.RUNNING
        handler.exitFile.text = ''

        then:
        !handler.checkIfCompleted()
        sleep 5_100
        handler.checkIfCompleted()
        handler.status == TaskHandler.Status.COMPLETED
        handler.task.exitStatus == Integer.MAX_VALUE

    }


    def testCheckIfTerminateEmptyWithLatency() {

        def task = new TaskRun()
        task.workDirectory = Files.createTempDirectory('testHandler')
        def config = Mock(TaskConfig)
        def executor = Mock(AbstractGridExecutor)

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
