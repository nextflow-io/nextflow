/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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
        def task = Mock(TaskRun)
        task.getName() >> 'tes-task'
        task.getWorkDir() >> Paths.get(".")
        task.getConfig() >> new TaskConfig(memory: '2GB', cpus: 4, disk: '10GB')
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


}
