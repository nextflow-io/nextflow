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
 * Execute a task script by running it on the SGE/OGE cluster
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class SgeExecutor extends AbstractGridExecutor {


    /*
     * Prepare the 'qsub' cmdline. The following options are used
     * - wd: define the job working directory
     * - terse: output just the job id on the output stream
     * - j: redirect qsub stdout and stderr to the same file (join)
     * - sync: wait for the job completion
     * -V: export the current environment
     */
    protected List<String> getSubmitCommandLine(TaskRun task) {

        final result = new ArrayList<String>()

        result << 'qsub'
        result << '-wd' << task.workDirectory?.toString()
        result << '-N' << "nf-${task.processor.name}-${task.index}"
        result << '-o' << "/dev/null"
        result << '-j' << 'y'
        result << '-terse'
        result << '-V'

        // the requested queue name
        if( taskConfig.queue ) {
            result << '-q'  << taskConfig.queue
        }

        // max task duration
        if( taskConfig.maxDuration ) {
            result << '-l' << "h_rt=${taskConfig.maxDuration.format('HH:mm:ss')}"
        }

        // task max memory
        if( taskConfig.maxMemory ) {
            result << '-l' << "virtual_free=${taskConfig.maxMemory.toString().replaceAll(/[\sB]/,'')}"
        }

        // -- at the end append the command script wrapped file name
        if( taskConfig.clusterOptions ) {
            result.addAll( getClusterOptionsAsList() )
        }

        // -- last entry to 'script' file name
        result << JOB_SCRIPT_FILENAME

        return result
    }

    @Override
    void killTask(jobId) {
       new ProcessBuilder( "qdel -j $jobId".split(' ') ).start()
    }
}
