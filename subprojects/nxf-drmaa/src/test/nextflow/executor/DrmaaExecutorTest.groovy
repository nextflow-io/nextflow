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
import nextflow.processor.TaskId
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import nextflow.processor.TaskStatus
import nextflow.trace.TraceRecord
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
        //def config = new TaskConfig([:])
        def executor = [:] as DrmaaExecutor
        def task = new TaskRun(id: TaskId.of(100), name: 'Hello', workDir: workDir, script: 'echo hello', config: [:])
        task.processor = Mock(TaskProcessor)
        task.processor.getProcessEnvironment() >> [:]
        task.processor.getSession() >> new nextflow.Session()
        task.processor.getConfig() >> [:]

        when:
        def handler = executor.createTaskHandler(task)
        then:
        handler instanceof DrmaaTaskHandler
        handler.workDir == workDir
        handler.taskName == 'nf-Hello'

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

        def task = new TaskRun(id: TaskId.of(1), name: 'hello', workDir: workDir, config: [queue: 'short'])
        def executor = [ getDrmaaSession: { drmaa } ] as DrmaaExecutor

        def handler = Spy(DrmaaTaskHandler, constructorArgs: [task, executor])
        when:
        handler.submit()
        then:
        1 * handler.createBashWrapper() >> null
        handler.status == TaskStatus.SUBMITTED
        handler.jobId == '12345'

        template.getWorkingDirectory() == workDir.toString()
        template.getRemoteCommand() == '/bin/bash'
        template.getArgs() == [handler.wrapperFile.toString()]
        template.getJoinFiles() == true
        template.getOutputPath() == ':' + workDir.resolve('.command.log')
        template.getNativeSpecification() == handler.getOptions()

        cleanup:
        workDir?.deleteDir()

    }

    def testHandlerGetOptions () {

        given:
        def workDir = Files.createTempDirectory('test')
        def executor = [:] as DrmaaExecutor
        def task = new TaskRun(id: TaskId.of(1), name: 'hello', workDir: workDir)
        def handler = new DrmaaTaskHandler(task, executor)
        def config = task.config = new TaskConfig()
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

        def expected = new TraceRecord()
        expected.task_id = TaskId.of(30)
        expected.native_id = '2000'
        expected.hash = '123abc'
        expected.name = 'hello (1)'
        expected.process = 'hello'
        expected.tag = 'seq1'
        expected.status = TaskStatus.SUBMITTED.toString()
        expected.exit = 99
        expected.submit = 1406264935000
        expected.start = 1406265009000
        expected.module = []
        expected.container = null
        expected.attempt = 1
        expected.workdir = workDir.toString()
        expected.script = null
        expected.scratch = null
        expected.queue = null
        expected.cpus = 1
        expected.time = null
        expected.memory = null
        expected.disk = null
        expected.env = null

        def task = [:] as TaskRun
        task.id = TaskId.of(30)
        task.workDir = workDir
        task.exitStatus = 99
        task.config = new TaskConfig(tag: 'seq1')
        task.processor = Mock(TaskProcessor)
        task.processor.getName() >> 'hello'
        task.processor.getSession() >> new nextflow.Session()
        task.processor.getProcessEnvironment() >> [:]
        task.metaClass.getName = { 'hello (1)' }
        task.metaClass.getHashLog =  {'123abc'}

        when:
        def handler = new DrmaaTaskHandler(task, Mock(DrmaaExecutor))
        handler.jobId = '2000'
        handler.status = TaskStatus.SUBMITTED
        handler.submitTimeMillis = 1406264935000
        handler.startTimeMillis = 1406265009000
        then:
        handler.getTraceRecord() == expected

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

