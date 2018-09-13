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

package nextflow.ga4gh.tes.executor

import java.nio.file.Paths

import nextflow.processor.TaskConfig
import nextflow.processor.TaskRun
import spock.lang.Specification
/**
 * @author Emilio Palumbo <emiliopalumbo@gmail.com>
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
