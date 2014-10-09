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

import nextflow.processor.ParallelTaskProcessor
import nextflow.processor.TaskConfig
import nextflow.processor.TaskHandler
import nextflow.processor.TaskRun
import nextflow.script.BaseScript
import nextflow.trace.TraceRecord
import org.ggf.drmaa.JobInfo
import org.ggf.drmaa.JobTemplate
import org.ggf.drmaa.Session
import spock.lang.Specification
import test.TestHelper
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DrmaaExecutorTest extends Specification {


    def testExecutorGetHandler() {

        given:
        def workDir = Files.createTempDirectory('temp')
        def config = new TaskConfig([:])
        def processor = [ getTaskConfig: {config}, getProcessEnvironment: {[:]} ] as ParallelTaskProcessor
        def executor = [:] as DrmaaExecutor
        def task = new TaskRun(id: 100, name: 'Hello', workDir: workDir, script: 'echo hello', processor: processor)

        when:
        def handler = executor.createTaskHandler(task)
        then:
        handler instanceof DrmaaTaskHandler
        handler.workDir == workDir.toFile()
        handler.taskName == 'nf-Hello'
        handler.wrapperFile.exists()

        cleanup:
        workDir?.deleteDir()

    }

    def testHandlerSubmit() {

        given:
        def workDir = Files.createTempDirectory('test')
        def template = TestHelper.proxyFor(JobTemplate)

        Session drmaa = Stub()
        drmaa.createJobTemplate() >> { template }
        drmaa.runJob(_) >> '12345'

        def task = new TaskRun(id:1, name: 'hello', workDir: workDir)
        def config = new TaskConfig([queue: 'short'])
        def executor = [ getDrmaaSession: { drmaa } ] as DrmaaExecutor

        def handler = new DrmaaTaskHandler(task, config, executor)
        when:
        handler.submit()

        then:
        handler.status == TaskHandler.Status.SUBMITTED
        handler.jobId == '12345'

        template.getWorkingDirectory() == workDir.toString()
        template.getRemoteCommand() == '/bin/bash'
        template.getArgs() == [handler.wrapperFile.toString()]
        template.getJoinFiles() == true
        template.getOutputPath() == ':/dev/null'
        template.getNativeSpecification() == handler.getOptions()

        cleanup:
        workDir?.deleteDir()

    }

    def testHandlerGetOptions () {

        given:
        def workDir = Files.createTempDirectory('test')
        def executor = [:] as DrmaaExecutor
        def task = new TaskRun(id:1, name: 'hello', workDir: workDir)
        def base = Mock(BaseScript)
        def config = new TaskConfig(base)
        def handler = new DrmaaTaskHandler(task, config, executor)

        when:
        config.queue = test_queue
        config.clusterOptions = test_opts
        config.time = test_time
        config.memory = test_mem
        config.cpus = test_cpus
        config.penv = test_penv

        then:
        handler.getOptions() == expected

        cleanup:
        workDir?.deleteDir()

        where:
        test_queue  | test_opts | test_mem  | test_time | test_cpus | test_penv || expected
        null        | null      | null      | null      | null      | null      || '-notify -b y'
        'short1'    | null      | null      | null      | null      | null      || '-notify -q short1 -b y'
        'long2'     | '-z abc'  | null      | null      | null      | null      || '-notify -q long2 -z abc -b y'
        'long2'     | '-z abc'  | '1M'      | null      | null      | null      || '-notify -q long2 -l virtual_free=1M -z abc -b y'
        'long2'     | '-z abc'  | '1G'      | '2h'      | null      | null      || '-notify -q long2 -l h_rt=02:00:00 -l virtual_free=1G -z abc -b y'
        'alpha'     | '-z abc'  | '2G'      | '3h'      | 2         | null      || '-notify -q alpha -l slots=2 -l h_rt=03:00:00 -l virtual_free=2G -z abc -b y'
        'delta'     | '-z abc'  | '2G'      | '3h'      | 2         | 'mpi'     || '-notify -q delta -pe mpi 2 -l h_rt=03:00:00 -l virtual_free=2G -z abc -b y'

    }

    def testHandlerGetTrace() {

        given:
        def workDir = Files.createTempDirectory('test')
        def usage = [
                start_time: '1406265009.0000',
                submission_time: '1406264935.0000',
                end_time: '1406265009.0000',
                maxvmem: '100.0000',
                cpu: '0.0040',
                io: '200.000'
        ]
        def jobInfo = [getResourceUsage:{ usage }] as JobInfo
        Session drmaa = Stub()
        drmaa.wait(_,0) >> jobInfo
        def executor = [ getDrmaaSession: { drmaa } ] as DrmaaExecutor
        def task = new TaskRun(id:30, name: 'hello', workDir: workDir, exitStatus: 99)
        def config = new TaskConfig([:])

        when:
        def handler = new DrmaaTaskHandler(task, config, executor)
        handler.jobId = '2000'
        handler.status = TaskHandler.Status.SUBMITTED
        then:
        handler.getTraceRecord() == new TraceRecord(
                taskId: 30,
                nativeId: '2000',
                name: 'hello',
                status: TaskHandler.Status.SUBMITTED,
                exit: 99,
                start: 1406265009000,
                submit: 1406264935000,
                complete: 1406265009000,
                cpu: '0.0040',
                mem: '100.0000')


        cleanup:
        workDir?.deleteDir()
    }

    def testToMillis() {

        expect:
        DrmaaTaskHandler.millis('xx') == 0
        DrmaaTaskHandler.millis(null) == 0

        //
        DrmaaTaskHandler.millis('1408691877.1200')    == 1408691877120
        DrmaaTaskHandler.millis('1409064132425.0000') == 1409064132425


    }

}

