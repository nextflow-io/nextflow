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

import nextflow.Session
import nextflow.processor.ParallelTaskProcessor
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
        def config = new TaskConfig()
        config.queue 'my-queue'
        config.maxMemory '2GB'
        config.maxDuration '3h'
        config.gridNativeOptions '-extra opt'
        config.name 'task'

        def src = new BaseScript() {
            @Override Object run() { return null  }
        }
        def executor = new SgeExecutor()
        def processor = new ParallelTaskProcessor(executor, new Session(), src, config, {})


        when:
        def task = new TaskRun(name: 'my-task', index: 9, processor: processor)
        task.workDirectory = new File('/abc')

        then:
        executor.getSubmitCommandLine(task) == 'qsub -wd /abc -N nf-task-9 -o /dev/null -j y -sync y -V -q my-queue -l h_rt=03:00:00 -l virtual_free=2G -extra opt .job.run'.split(' ') as List

    }

}
