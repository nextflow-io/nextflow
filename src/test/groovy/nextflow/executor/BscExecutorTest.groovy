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

import java.nio.file.Path
import java.nio.file.Paths

import nextflow.processor.TaskConfig
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class BscExecutorTest extends Specification {

    def testHeaders() {

        setup:
        // mock process
        def proc = Mock(TaskProcessor)
        def config = new TaskConfig()

        // LSF executor
        def executor = [:] as BscExecutor

        when:
        // process name
        proc.getName() >> 'task'
        // config
        config.queue = 'bsc_ls'
        config.clusterOptions = "-x 1 -R \"span[ptile=2]\""
        config.cpus = '2'
        config.time = '1h 30min'
        config.memory = '8GB'

        // task object
        def task = new TaskRun()
        task.config = config
        task.processor = proc
        task.workDir = Paths.get('/scratch')
        task.index = 1

        then:
        executor.getHeaders(task) == '''
                #BSUB -cwd /scratch
                #BSUB -o /dev/null
                #BSUB -q bsc_ls
                #BSUB -n 2
                #BSUB -R span[hosts=1]
                #BSUB -W 01:30
                #BSUB -M 4096
                #BSUB -J nf-task_1
                #BSUB -x 1
                #BSUB -R span[ptile=2]
                '''
                .stripIndent().leftTrim()

    }


    def testWrapString() {

        given:
        def executor = [:] as BscExecutor

        expect:
        executor.wrap('') == ''
        executor.wrap('hello') == 'hello'
        executor.wrap('span[ptile=number]') == '"span[ptile=number]"'
        executor.wrap('hello world') == '"hello world"'

        executor.pipeLauncherScript() == true
        executor.getSubmitCommandLine(Mock(TaskRun), Mock(Path)) == ['bsub']
    }

}
