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

/**
 * Processor for LSF resource manager (DRAFT)
 *
 * See http://en.wikipedia.org/wiki/Platform_LSF
 *
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@InheritConstructors
class LsfTaskProcessor extends GenericGridProcessor {


    protected String bsubCmdLine

    /**
     * Extra options appended to the generated 'qsub' command line
     */
    LsfTaskProcessor bsubCmdLine( String cmdLine ) {
        this.bsubCmdLine = cmdLine
        return this
    }

    @Override
    protected List<String> getSubmitCommandLine(TaskRun task) {

        final result = new ArrayList<String>()

        result << 'bsub'
        result << '-cwd' << task.workDirectory
        result << '-K'              // sync mode i.e. wait for termination before exit
        result << '-V'

        // add other parameters (if any)
        if(queue) {
            result << '-q'  << queue
        }

        if( maxDuration ) {
            result << '-l' << "h_rt=${maxDuration.format('HH:mm:ss')}"
        }

        if( maxMemory ) {
            result << '-l' << "virtual_free=${maxMemory.toString().replaceAll(/[\sB]/,'')}"
        }

        // -- the job name
        result << '-J' << "nf-${name}-${task.index}"

        // -- at the end append the command script wrapped file name
        if ( bsubCmdLine ) {
            if( bsubCmdLine instanceof Collection ) {
                result.addAll( bsubCmdLine as Collection )
            }
            else {
                result.addAll( bsubCmdLine.toString().split(' ') )
            }
        }

        // -- last entry to 'script' file name
        result << '<' << JOB_SCRIPT_FILENAME

        return result

    }
}
