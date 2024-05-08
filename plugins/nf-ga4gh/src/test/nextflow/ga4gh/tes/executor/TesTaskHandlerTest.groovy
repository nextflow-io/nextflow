/*
 * Copyright 2013-2024, Seqera Labs
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

package nextflow.ga4gh.tes.executor

import nextflow.ga4gh.tes.client.api.TaskServiceApi
import nextflow.ga4gh.tes.client.model.TesCreateTaskResponse
import nextflow.ga4gh.tes.client.model.TesTask
import nextflow.processor.TaskStatus

import java.nio.file.Path
import java.nio.file.Paths

import nextflow.processor.TaskConfig
import nextflow.processor.TaskRun
import spock.lang.Specification
/**
 * @author Emilio Palumbo <emilio.palumbo@crg.eu>
 */
class TesTaskHandlerTest extends Specification {

    def 'should set resources' () {

        given:
        def executor = Mock(TesExecutor)
        def task = Mock(TaskRun) {
            getInputFilesMap() >> [:]
            getOutputFilesNames() >> []
            getOutputsByType(_) >> [:]
        }
        task.getName() >> 'tes-task'
        task.getWorkDir() >> Paths.get(".")
        task.getConfig() >> new TaskConfig(memory: '2GB', cpus: 4, disk: '10GB')
        task.getContainer() >> 'foo'
        def handler = new TesTaskHandler(task, executor)


        when:
        def t = handler.newTesTask()

        then:
        t.getResources().cpuCores == 4
        t.getResources().ramGb == 2
        t.getResources().diskGb == 10
        t.getResources().preemptible == null
        t.getResources().zones == null

    }

    def 'should submit job' () {

        given:
        def executor = Mock(TesExecutor)
        def task = Mock(TaskRun)
        def bashBuilder = Mock(TesBashBuilder)
        def client = Mock(TaskServiceApi)
        def handler = Spy(TesTaskHandler)
        handler.@client = client
        handler.@executor = executor
        handler.task = task

        def req = Mock(TesTask)
        def resp = Mock(TesCreateTaskResponse)

        when:
        handler.submit()
        then:
        1 * executor.getRemoteBinDir() >> Path.of("/work/bin")
        1 * handler.newTesBashBuilder(task, "/work/bin") >> bashBuilder
        1 * bashBuilder.build() >> null
        1 * handler.newTesTask() >> req
        1 * client.createTask(req) >> resp
        1 * resp.getId() >> '12345'

        handler.status == TaskStatus.SUBMITTED
        handler.requestId == '12345'
    }

}
