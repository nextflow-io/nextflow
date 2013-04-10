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

package nextflow.processor
import groovy.transform.InheritConstructors
import groovy.util.logging.Slf4j
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@InheritConstructors
class NopeTaskProcessor extends AbstractTaskProcessor {

    @Override
    protected void launchTask( TaskDef task ) {

        task.workDirectory = new File('.').absoluteFile
        task.status = TaskDef.Status.TERMINATED
        task.exitCode = 0
        task.output = task.script   // return the script itself as output

    }

    @Override
    protected collectResultFile(TaskDef task, String name) {
        return new File(task.workDirectory, name)
    }
}
