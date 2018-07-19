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

import java.nio.file.Paths

import nextflow.processor.TaskConfig
import nextflow.processor.TaskRun
import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class SlurmTaskHandlerTest extends Specification {


    def 'should append --cluster to the list' () {
        given:
        def exec  = Mock(SlurmExecutor)
        def task = Mock(TaskRun)
        SlurmTaskHandler handler

        when:
        handler = new SlurmTaskHandler(task, exec)
        then:
        task.workDir >> Paths.get('/foo')
        task.getConfig() >> [:]
        handler.getCluster() == null

        when:
        handler = new SlurmTaskHandler(task, exec)
        then:
        task.workDir >> Paths.get('/foo')
        task.getConfig() >> [clusterOptions:'--cluster foo']
        handler.getCluster() == 'foo'
    }


    @Unroll
    def 'should fetch cluster option with string=#STR' () {

        given:
        def h = new SlurmTaskHandler()

        expect:
        h.fetchClusterOption(STR) == EXPECTED

        where:
        STR                    | EXPECTED
        null                   | null
        ''                     | null
        '--this --that '       | null          
        '--cluster='           | null
        '--cluster foo'        | 'foo'
        '--cluster=bar --that' | 'bar'

    }

    def 'should return cluster options' () {

        given:
        def task = Mock(TaskRun)
        def exec = Mock(SlurmExecutor)
        def conf = Mock(TaskConfig)

        when:
        def opts = new SlurmTaskHandler().getQueueOpts()
        then:
        opts == Collections.emptyMap()

        when:
        opts = new SlurmTaskHandler(task,exec).getQueueOpts()
        then:
        task.getWorkDir() >> Paths.get('/some/dir')
        task.getConfig() >> conf
        1 * conf.getProperty('queue') >> 'foo'
        1 * conf.getProperty('clusterOptions') >> '--cluster bar'

        opts == [queue:'foo', cluster:'bar']
    }

}
