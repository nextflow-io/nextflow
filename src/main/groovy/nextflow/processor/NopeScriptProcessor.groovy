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
class NopeScriptProcessor extends AbstractScriptProcessor {

    @Override
    protected void runScript( def script, TaskDef task ) {
        log.debug "Running script: $script"

        task.workDirectory = new File('.').absoluteFile
        task.status = TaskDef.Status.TERMINATED
        task.exitCode = 0
        task.output = script ?. toString()  // return the script itself as output

    }

    @Override
    protected List<File> collectResultFile(File scratchPath, String name) {
        return [ new File(scratchPath, name) ]
    }
}
