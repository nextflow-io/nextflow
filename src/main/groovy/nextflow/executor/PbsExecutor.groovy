/*
 * Copyright (c) 2013-2015, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2015, Paolo Di Tommaso and the respective authors.
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
import groovy.util.logging.Slf4j
import nextflow.processor.TaskRun
/**
 * Implements a executor for PBS/Torque cluster
 *
 * See http://www.pbsworks.com
 */
@Slf4j
class PbsExecutor extends AbstractGridExecutor {

    /**
     * Gets the directives to submit the specified task to the cluster for execution
     *
     * @param task A {@link TaskRun} to be submitted
     * @param result The {@link List} instance to which add the job directives
     * @return A {@link List} containing all directive tokens and values.
     */
    protected List<String> getDirectives( TaskRun task, List<String> result ) {
        assert result !=null

        result << '-d' << task.workDir.toString()
        result << '-N' << getJobNameFor(task)
        result << '-o' << task.workDir.resolve(TaskRun.CMD_LOG).toString()
        result << '-j' << 'oe'
        result << '-V' << ''

        // the requested queue name
        if( task.config.queue ) {
            result << '-q'  << (String)task.config.queue
        }

        if( task.config.cpus > 1 ) {
            result << '-l' << "nodes=1:ppn=${task.config.cpus}"
        }

        // max task duration
        if( task.config.time ) {
            final duration = task.config.getTime()
            result << "-l" << "walltime=${duration.format('HH:mm:ss')}"
        }

        // task max memory
        if( task.config.memory ) {
            result << "-l" << "mem=${task.config.memory.toString().replaceAll(/[\s]/,'').toLowerCase()}"
        }

        // -- at the end append the command script wrapped file name
        if( task.config.clusterOptions ) {
            result << task.config.clusterOptions.toString() << ''
        }

        return result
    }

    /**
     * The command line to submit this job
     *
     * @param task The {@link TaskRun} instance to submit for execution to the cluster
     * @param scriptFile The file containing the job launcher script
     * @return A list representing the submit command line
     */
    List<String> getSubmitCommandLine(TaskRun task, Path scriptFile ) {
        [ 'qsub', scriptFile.getName() ]
    }

    protected String getHeaderToken() { '#PBS' }

    /**
     * Parse the string returned by the {@code qsub} command and extract the job ID string
     *
     * @param text The string returned when submitting the job
     * @return The actual job ID string
     */
    @Override
    def parseJobId( String text ) {
        // return always the last line
        def result = text?.trim()
        if( result ) {
            return result.tokenize('.').get(0)
        }

        throw new IllegalStateException("Invalid PBS/Torque submit response:\n$text\n\n")
    }


    @PackageScope
    String getKillCommand() { 'qdel' }

    @Override
    protected List<String> queueStatusCommand(Object queue) {
        def result = ['qstat']
        if( queue )
            result << queue.toString()

        return result
    }

    static private Map DECODE_STATUS = [
            'C': QueueStatus.DONE,
            'R': QueueStatus.RUNNING,
            'Q': QueueStatus.PENDING,
            'H': QueueStatus.HOLD,
            'S': QueueStatus.HOLD
    ]

    @Override
    protected Map<?, QueueStatus> parseQueueStatus(String text) {

        def result = [:]
        text?.eachLine{ String row, int index ->
            if( index< 2 ) return
            def cols = row.split(/\s+/)
            if( cols.size()>5 ) {
                result.put( parseJobId(cols[0]), DECODE_STATUS[cols[4]] ?: AbstractGridExecutor.QueueStatus.UNKNOWN )
            }
        }

        return result
    }


}
