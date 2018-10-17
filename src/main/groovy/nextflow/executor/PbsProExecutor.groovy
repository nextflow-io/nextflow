/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
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
import nextflow.processor.TaskRun



/**
 * Implements a executor for PBSPro cluster
 * Tested with version:
 * - 14.2.4
 * - 19.0.0
 * See http://www.pbspro.org
 *
 * @author Lorenz Gerber <lorenzottogerber@gmail.com>
 */
@Slf4j
class PbsProExecutor extends PbsExecutor {

    /**
     * Gets the directives to submit the specified task to the cluster for execution
     *
     * @param task A {@link TaskRun} to be submitted
     * @param result The {@link List} instance to which add the job directives
     * @return A {@link List} containing all directive tokens and values.
     */
    @Override
    protected List<String> getDirectives(TaskRun task, List<String> result ) {
        assert result !=null

        result << '-N' << getJobNameFor(task)
        result << '-o' << quote(task.workDir.resolve(TaskRun.CMD_LOG))
        result << '-j' << 'oe'

        // the requested queue name
        if( task.config.queue ) {
            result << '-q'  << (String)task.config.queue
        }

        if( task.config.cpus == 1 && task.config.memory ) {
            result << '-l' << "select=mem=${task.config.memory.toString().replaceAll(/[\s]/,'').toLowerCase()}"
        }

        if( task.config.cpus > 1 && !task.config.memory ) {
            result << '-l' << "select=1:ncpus=${task.config.cpus}"
        }


        if( task.config.cpus > 1 && task.config.memory ) {
            result << '-l' << "select=1:ncpus=${task.config.cpus}:mem=${task.config.memory.toString().replaceAll(/[\s]/,'').toLowerCase()}"
        }

        // max task duration
        if( task.config.time ) {
            final duration = task.config.getTime()
            result << "-l" << "walltime=${duration.format('HH:mm:ss')}"
        }

        // -- at the end append the command script wrapped file name
        if( task.config.clusterOptions ) {
            result << task.config.clusterOptions.toString() << ''
        }

        return result
    }

    @Override
    protected List<String> queueStatusCommand(Object queue) {
        String cmd = 'qstat -f'
        if( queue ) cmd += ' ' + queue
        return ['sh','-c', "$cmd | egrep '(Job Id:|job_state =)'".toString()]
    }

    static private Map DECODE_STATUS = [
            'C': QueueStatus.DONE,
            'R': QueueStatus.RUNNING,
            'Q': QueueStatus.PENDING,
            'H': QueueStatus.HOLD,
            'S': QueueStatus.HOLD
    ]

}
