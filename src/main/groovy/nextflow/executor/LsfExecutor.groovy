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

import nextflow.processor.TaskRun

/**
 * Processor for LSF resource manager (DRAFT)
 *
 * See http://en.wikipedia.org/wiki/Platform_LSF
 *
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class LsfExecutor extends AbstractGridExecutor {

    @Override
    protected List<String> getSubmitCommandLine(TaskRun task) {

        final result = new ArrayList<String>()

        result << 'bsub'
        result << '-K'    // sync mode i.e. wait for termination before exit
        result << '-cwd' << task.workDirectory?.toString()
        result << '-o' << JOB_OUT_FILENAME

        // add other parameters (if any)
        if( taskConfig.queue ) {
            result << '-q'  << taskConfig.queue
        }

        // -- the job name
        result << '-J' << "nf-${task.processor.name}-${task.index}"

        // -- at the end append the command script wrapped file name
        if( taskConfig.gridNativeOptions ) {
            result.addAll( getGridNativeOptionsAsList() )
        }

        // -- last entry to 'script' file name
        result << "./$JOB_SCRIPT_FILENAME"

        return result

    }

    def submitJob( TaskRun task, File runnerFile, File cmdOutFile ) {
        // note: LSF requires the job script file to be executable
        runnerFile.setExecutable(true)

        // now invoke the default method
        super.submitJob(task, runnerFile, cmdOutFile)
    }
}
