/*
 * Copyright 2019, Genome Research Limited
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

package nextflow.wr.executor

import java.nio.file.Path
import java.nio.file.Paths
import java.nio.file.Files

import nextflow.Session
import nextflow.processor.TaskConfig
import nextflow.processor.TaskRun
import nextflow.processor.TaskProcessor
import nextflow.processor.BatchContext
import nextflow.wr.client.WrRestApi
import spock.lang.Specification
/**
 * @author Sendu Bala <sb10@sanger.ac.uk>
 * Based on TesTaskHandlerTest Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class WrTaskHandlerTest extends Specification {

    def 'should return job ids as a string' () {
        given:
        def executor = Mock(WrExecutor)
        def task = Mock(TaskRun)
        def handler = new WrTaskHandler(task, executor)

        when:
        Collection ids = ['1','2','3']
        def str = handler.jobIdsToString(ids)

        then:
        str == '1,2,3'

        when:
        ids = ['1','2','3','4','5','6','7','8','9','10','11']
        str = handler.jobIdsToString(ids)

        then:
        str == '1,2,3,4,5,6,7,8,9,10, ... other 1 omitted'
    }

    def 'should return job details unbatched' () {
        given:
        def executor = Mock(WrExecutor)
        def task = Mock(TaskRun)
        Map expectedJob = ["Key": "one"]
        List<Map> jobs = [expectedJob]
        WrRestApi client = Mock {
            1 * status("one") >> jobs
        }
        def handler = new WrTaskHandler(task, executor)
        handler.client = client

        when:
        def job = handler.getJob("one")

        then:
        job == expectedJob
    }

    def 'should return job details batched' () {
        given:
        def executor = Mock(WrExecutor)
        def task = Mock(TaskRun)
        Map expectedJob1 = ["Key": "one"]
        Map expectedJob2 = ["Key": "two"]
        List<Map> jobs = [expectedJob1, expectedJob2]
        WrRestApi client = Mock {
            1 * status("one,two") >> jobs
        }
        def handler1 = new WrTaskHandler(task, executor)
        handler1.client = client
        handler1.submitted("one")
        def handler2 = new WrTaskHandler(task, executor)
        handler2.client = client
        handler2.submitted("two")

        BatchContext c = Spy(BatchContext)

        when:
        handler1.batch(c)
        handler2.batch(c)
        def job1 = handler1.getJob("one")
        def job2 = handler2.getJob("two")

        then:
        1 * c.collect('one')
        1 * c.collect("two")
        1 * c.getBatchFor('one', 1000)
        1 * c.get('two')
        job1 == expectedJob1
        job2 == expectedJob2
    }

    def 'should return if running or completed' () {
        given:
        def executor = Mock(WrExecutor)
        def task = [name: 'Hello 1'] as TaskRun
        List<Map> jobs = [
            ["Key": "one", "State": "running"],
            ["Key": "two", "State": "complete", "Exited": true, "Exitcode": "3"]
        ]
        WrRestApi client = Mock {
            4 * status(_) >> jobs
        }
        def handler = new WrTaskHandler(task, executor)
        handler.client = client

        when:
        def running = handler.checkIfRunning()
        def complete = handler.checkIfCompleted()

        then:
        assert !running
        assert !complete

        when:
        handler.submitted("one")
        running = handler.checkIfRunning()
        complete = handler.checkIfCompleted()

        then:
        assert running
        assert !complete

        when:
        handler.submitted("two")
        running = handler.checkIfRunning()
        complete = handler.checkIfCompleted()

        then:
        assert running
        assert complete
        task.exitStatus == 3
    }

    def 'should be able to kill a job' () {
        given:
        def executor = Mock(WrExecutor)
        def task = Mock(TaskRun)
        def client = Mock(WrRestApi)
        def handler = new WrTaskHandler(task, executor)
        handler.client = client
        String jobId = "id"

        when:
        handler.submitted(jobId)
        handler.kill()

        then:
        1 * client.cancel(jobId)
    }

    def 'should be able to get submit args' () {
        given:
        def folder = Files.createTempDirectory('test')
        def executor = Mock(WrExecutor)
        def task = Spy([
            name: 'Hello 1',
            workDir: folder,
            script: 'echo Hello world!'
        ] as TaskRun)
        def config = Mock(TaskConfig)
        task.config = config
        def processor = Mock(TaskProcessor)
        def session = Mock(Session)
        def client = Mock(WrRestApi)
        def handler = new WrTaskHandler(task, executor)
        handler.client = client

        when:
        def returned = handler.submitArgs()

        then:
        task.getEnvironment() >> null
        task.getCondaEnv() >> null
        task.getContainerConfig() >> null
        task.isContainerNative() >> false
        task.getProcessor() >> processor
        processor.getSession() >> session
        processor.getConfig() >> config
        session.getWorkDir() >> folder
        returned.size() == 3
        returned[0] == "/bin/bash $folder/.command.run"
        returned[1] == task
        assert returned[2] instanceof WrFileCopyStrategy 

        cleanup:
        folder.deleteDir()
    }

}
