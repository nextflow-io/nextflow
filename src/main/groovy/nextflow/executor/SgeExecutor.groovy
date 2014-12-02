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

import groovy.transform.PackageScope
import nextflow.processor.TaskRun
import nextflow.util.Duration

/**
 * Execute a task script by running it on the SGE/OGE cluster
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class SgeExecutor extends AbstractGridExecutor {

    /**
     * Gets the directives to submit the specified task to the cluster for execution
     *
     * @param task A {@link TaskRun} to be submitted
     * @param result The {@link List} instance to which add the job directives
     * @return A {@link List} containing all directive tokens and values.
     */
    protected List<String> getDirectives(TaskRun task, List<String> result) {

        result << '-wd' << task.workDir?.toString()
        result << '-N' << getJobNameFor(task)
        result << '-o' << '/dev/null'
        result << '-j' << 'y'
        result << '-terse'
        result << '-V'
        /*
         * By using command line option -notify SIGUSR1 will be sent to your script prior to SIGSTOP
         * and SIGUSR2 will be sent to your script prior to SIGKILL
         */
        result << '-notify'

        // the requested queue name
        if( taskConfig.queue ) {
            result << '-q' << (taskConfig.queue as String)
        }

        //number of cpus for multiprocessing/multi-threading
        if( taskConfig.cpus ) {
            if ( taskConfig.penv ) {
                result << "-pe" << taskConfig.penv.toString() << taskConfig.cpus.toString()
            } else {
                result << "-l" << "slots=${taskConfig.cpus}"
            }
        }

        // max task duration
        if( taskConfig.time ) {
            final time = taskConfig.time as Duration
            result << "-l" << "h_rt=${time.format('HH:mm:ss')}"
        }

        // task max memory
        if( taskConfig.memory ) {
            result << "-l" << "virtual_free=${taskConfig.memory.toString().replaceAll(/[\sB]/,'')}"
        }

        // -- at the end append the command script wrapped file name
        if( taskConfig.clusterOptions ) {
            result.addAll( getClusterOptionsAsList() )
        }

        return result
    }


    /*
     * Prepare the 'qsub' cmdline. The following options are used
     * - wd: define the job working directory
     * - terse: output just the job id on the output stream
     * - j: redirect qsub stdout and stderr to the same file (join)
     * - sync: wait for the job completion
     * -V: export the current environment
     */
    List<String> getSubmitCommandLine(TaskRun task, Path scriptFile ) {

        final result = ['qsub']

        // -- adds the jobs directives to the command line
        getDirectives(task, result)

        // -- last entry to 'script' file name
        result << scriptFile.getName()

        return result
    }

    /**
     * Parse the string returned by the {@code qsub} command and extract the job ID string
     *
     * @param text The string returned when submitting the job
     * @return The actual job ID string
     */
    @Override
    def parseJobId( String text ) {
        // return always the last line
        def lines = text.trim().readLines()
        def id = lines[-1].trim()
        if( id?.toString()?.isInteger() )
            return id

        throw new IllegalStateException("Invalid SGE submit response:\n$text\n\n")
    }


    @PackageScope
    List<String> killTaskCommand(jobId) {
        ['qdel', '-j', jobId?.toString()]
    }

    @Override
    protected List<String> queueStatusCommand(Object queue) {
        def result = ['qstat']
        if( queue )
            result << '-q' << queue.toString()

        return result
    }

    static protected Map DECODE_STATUS = [
            'r': QueueStatus.RUNNING,
            'qw': QueueStatus.PENDING,
            'hqw': QueueStatus.HOLD,
            'Eqw': QueueStatus.ERROR
    ]

    @Override
    protected Map<?, QueueStatus> parseQueueStatus(String text) {

        def result = [:]
        text?.eachLine{ String row, int index ->
            if( index< 2 ) return
            def cols = row.trim().split(/\s+/)
            if( cols.size()>5 ) {
                result.put( cols[0], DECODE_STATUS[cols[4]] )
            }
        }

        return result
    }

}
