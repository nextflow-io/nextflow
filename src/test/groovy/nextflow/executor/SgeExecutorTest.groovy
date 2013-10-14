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
import java.nio.file.Paths

import nextflow.processor.TaskConfig
import nextflow.processor.TaskRun
import nextflow.script.BaseScript
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class SgeExecutorTest extends Specification {

    def 'test qsub cmd line' () {

        setup:
        def script = Mock(BaseScript)
        def config = new TaskConfig(script)
        config.queue 'my-queue'
        config.maxMemory '2GB'
        config.maxDuration '3h'
        config.clusterOptions '-extra opt'
        config.name 'task'


        def executor = [:] as SgeExecutor
        executor.taskConfig = config

        when:
        def wrapper = Paths.get('.job.sh')
        def task = new TaskRun(name: 'task x')
        task.workDirectory = Paths.get('/abc')

        then:
        executor.getSubmitCommandLine(task,wrapper) == 'qsub -wd /abc -N nf-task_x -o /dev/null -j y -terse -V -q my-queue -l h_rt=03:00:00 -l virtual_free=2G -extra opt .job.sh'.split(' ') as List

    }


    def testParseJobId() {

        when:
        def executor = [:] as SgeExecutor
        def textToParse = '''
            blah blah ..
            .. blah blah
            6472
            '''
        then:
        executor.parseJobId(textToParse) == '6472'
    }

    def testKillTaskCommand() {

        when:
        def executor = [:] as SgeExecutor
        then:
        executor.killTaskCommand(123) == ['qdel', '-j', '123'] as String[]

    }

}
