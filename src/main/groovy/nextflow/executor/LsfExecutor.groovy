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

import groovy.transform.InheritConstructors
import nextflow.processor.TaskRun

/**
 * Processor for LSF resource manager (DRAFT)
 *
 * See http://en.wikipedia.org/wiki/Platform_LSF
 *
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@InheritConstructors
class LsfExecutor extends GenericGridExecutor {

    @Override
    protected List<String> getSubmitCommandLine(TaskRun task) {

        final result = new ArrayList<String>()

        result << 'bsub'
        result << '-cwd' << task.workDirectory?.toString()
        result << '-K'              // sync mode i.e. wait for termination before exit
        result << '-V'

        // add other parameters (if any)
        if( taskConfig.queue ) {
            result << '-q'  << taskConfig.queue
        }

        if( taskConfig.maxDuration ) {
            result << '-l' << "h_rt=${taskConfig.maxDuration.format('HH:mm:ss')}"
        }

        if( taskConfig.maxMemory ) {
            result << '-l' << "virtual_free=${taskConfig.maxMemory.toString().replaceAll(/[\sB]/,'')}"
        }

        // -- the job name
        result << '-J' << "nf-${task.processor.name}-${task.index}"

        // -- at the end append the command script wrapped file name
        if ( taskConfig.nativeGridOptions ) {
            if( taskConfig.nativeGridOptions instanceof Collection ) {
                result.addAll( taskConfig.nativeGridOptions as Collection )
            }
            else {
                result.addAll( taskConfig.nativeGridOptions.toString().split(' ') )
            }
        }

        // -- last entry to 'script' file name
        result << '<' << JOB_SCRIPT_FILENAME

        return result

    }
}
