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
import groovy.util.logging.Slf4j
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun

/**
 * Dummy executor, only for test purpose
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class NopeExecutor extends ExecutionStrategy {

    @Override
    void launchTask( TaskProcessor processor, TaskRun task ) {

        task.workDirectory = new File('.').absoluteFile
        task.status = TaskRun.Status.TERMINATED
        task.exitCode = 0
        task.output = task.script   // return the script itself as output

    }

    @Override
    def getStdOutFile(TaskRun task) {
        return task.script
    }

}
