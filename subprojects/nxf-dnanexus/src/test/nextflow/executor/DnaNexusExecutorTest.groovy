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

import nextflow.fs.dx.api.DxApi
import nextflow.processor.TaskConfig
import nextflow.processor.TaskHandler
import nextflow.processor.TaskRun
import nextflow.script.BaseScript
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DnaNexusExecutorTest extends Specification {

    def testCreateInputObject() {

        given:
        def exec = [:] as DnaNexusExecutor
        when:
        def obj  = exec.createInputObject( [a:1, b:2], 'dx_cc2.8xlarge' )
        then:
        obj.input == [a:1, b:2]
        obj.function == 'process'
        obj.systemRequirements.process.instanceType == 'dx_cc2.8xlarge'

    }


    def testHandlerSubmit() {

        given:
        def api = Mock(DxApi)
        def folder = Files.createTempDirectory('tst-path')
        def task = Mock(TaskRun)
        task.getWorkDir() >> folder
        def exec = Mock(DnaNexusExecutor)
        def script = Mock(BaseScript)

        def params = [file1:'abc', file2: 'xxx']
        // note: the *instanceType* configured must used in the submit method
        def config = new TaskConfig(script)
        config.instanceType = 'dx_m1.super'
        // define the handler
        def handler = new DxTaskHandler(task, config, exec, params, api);

        // constraints
        and:
        1 * exec.createInputObject(params,'dx_m1.super') >> [function:'process']
        1 * api.jobNew(_) >> "job-xyz"

        when:
        handler.submit()
        then:
        handler.processJobId == "job-xyz"

        cleanup:
        folder?.deleteDir()
    }


    def testKill() {

        given:
        def api = Mock(DxApi)

        def folder = Files.createTempDirectory('tst-path')
        def task = Mock(TaskRun)
        task.getWorkDir() >> folder

        def exec = Mock(DnaNexusExecutor)
        def config = Mock(TaskConfig)
        def handler = new DxTaskHandler(task, config, exec, null, api);
        handler.processJobId = 'job-123'

        and:
        1 * api.jobTerminate('job-123')

        when:
        handler.kill()
        then:
        noExceptionThrown()

        cleanup:
        folder?.deleteDir()
    }


    def testCheckStatus() {

        given:
        def api = Mock(DxApi)

        def folder = Files.createTempDirectory('tst-path')
        def task = Mock(TaskRun)
        task.getWorkDir() >> folder

        def exec = Mock(DnaNexusExecutor)
        def config = Mock(TaskConfig)
        def handler = new DxTaskHandler(task, config, exec, null, api);
        handler.processJobId = 'job-312'

        and:
        1 * api.jobDescribe('job-312') >> [alpha:1, beta:2]

        when:
        handler.checkStatus() ==  [alpha:1, beta:2]
        then:
        noExceptionThrown()

        cleanup:
        folder?.deleteDir()
    }


    def testHandleCheckIfRunning() {

        setup:
        def api = Mock(DxApi)

        def folder = Files.createTempDirectory('tst-path')
        def task = Mock(TaskRun)
        task.getWorkDir() >> folder

        def config = Mock(TaskConfig)
        def exec = Mock(DnaNexusExecutor)
        def handler = new DxTaskHandler(task, config, exec, [:], api);
        handler.metaClass.checkStatus = { return [state:'runnable'] }
        when:
        handler.processJobId = '123'
        handler.status = TaskHandler.Status.SUBMITTED
        then:
        handler.checkIfRunning()

        cleanup:
        folder?.deleteDir()
    }

    def testTaskHandlerCheckIfTerminated() {

        setup:
        def api = Mock(DxApi)

        def work = Files.createTempDirectory('testWorkDir')
        def task = Mock(TaskRun)
        task.getWorkDir() >> work

        def config = Mock(TaskConfig)
        def exec = Mock(DnaNexusExecutor)

        when:
        def handler = new DxTaskHandler(task, config, exec, [:], api);
        handler.metaClass.checkStatus = { return [state:'running'] }
        handler.processJobId = '123'
        handler.status = TaskHandler.Status.RUNNING
        then:
        !handler.checkIfCompleted()
        handler.status == TaskHandler.Status.RUNNING

        when:
        def task2 = new TaskRun()
        task2.workDir = Files.createTempDirectory('testHandler')
        task2.workDir.resolve(TaskRun.CMD_OUTFILE).text = 'Task says Hola'
        handler = new DxTaskHandler(task2, config, exec, [:], api);
        handler.metaClass.checkStatus = { return [state:'done', output:[exit_code:33]] }
        handler.processJobId = '123'
        handler.status = TaskHandler.Status.RUNNING
        then:
        handler.checkIfCompleted()
        handler.status == TaskHandler.Status.COMPLETED
        task2.exitStatus == 33
        task2.stdout == 'Task says Hola'


        cleanup:
        work?.deleteDir()
        task2?.workDir?.deleteDir()
    }




}