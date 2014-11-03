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

import groovy.transform.InheritConstructors
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.processor.TaskRun
import nextflow.util.Duration


@Slf4j
@InheritConstructors
class PbsExecutor extends AbstractGridExecutor {

    List<String> getSubmitCommandLine(TaskRun task, Path scriptFile ) {

        final result = new ArrayList<String>()

        result << 'qsub'
        result << '-d' << task.workDir?.toString()
        result << '-N' << getJobNameFor(task)
        result << '-o' << '/dev/null'
        result << '-e' << '/dev/null'
        result << '-V'

        // the requested queue name
        if( taskConfig.queue ) {
            result << '-q'  << (String)taskConfig.queue
        }

        if( taskConfig.cpus ) {
            result << '-l' << "nodes=1:ppn=${taskConfig.cpus}"
        }

        // max task duration
        if( taskConfig.time ) {
            final duration = taskConfig.time as Duration
            result << "-l" << "walltime=${duration.format('HH:mm:ss')}"
        }

        // task max memory
        if( taskConfig.memory ) {
            result << "-l" << "mem=${taskConfig.memory.toString().replaceAll(/[\s]/,'').toLowerCase()}"
        }

        // -- at the end append the command script wrapped file name
        if( taskConfig.clusterOptions ) {
            result.addAll( getClusterOptionsAsList() )
        }

        // -- last entry to 'script' file name
        result << scriptFile.getName()

        return result
    }

    @Override
    def parseJobId( String text ) {
        // return always the last line
        def result = text?.trim()
        if( result )
            return result

        throw new IllegalStateException("Invalid PBS/Torque submit response:\n$text\n\n")
    }


    @PackageScope
    List<String> killTaskCommand(jobId) {
        ['qdel', jobId?.toString()]
    }

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
                result.put( cols[0], DECODE_STATUS[cols[4]] ?: AbstractGridExecutor.QueueStatus.UNKNOWN )
            }
        }

        return result
    }


}