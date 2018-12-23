/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package nextflow.executor

import groovy.util.logging.Slf4j
import nextflow.processor.TaskRun

/**
 * Implements a executor for PBSPro cluster executor
 *
 * Tested with version:
 * - 14.2.4
 * - 19.0.0
 * See http://www.pbspro.org
 *
 * @author Lorenz Gerber <lorenzottogerber@gmail.com>
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
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

        def res = []
        if( task.config.cpus > 1 ) {
            res << "ncpus=${task.config.cpus}".toString()
        }
        if( task.config.memory ) {
            res << "mem=${task.config.getMemory().getMega()}mb".toString()
        }
        if( res ) {
            result << '-l' << "select=1:${res.join(':')}".toString()
        }

        // max task duration
        if( task.config.time ) {
            final duration = task.config.getTime()
            result << "-l" << "walltime=${duration.format('HH:mm:ss')}".toString()
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
        return ['bash','-c', "set -o pipefail; $cmd | { egrep '(Job Id:|job_state =)' || true; }".toString()]
    }

    // see https://www.pbsworks.com/pdfs/PBSRefGuide18.2.pdf
    // table 8.1
    static private Map DECODE_STATUS = [
            'F': QueueStatus.DONE,      // job is finished
            'E': QueueStatus.RUNNING,   // job is exiting (therefore still running)
            'R': QueueStatus.RUNNING,   // job is running 
            'Q': QueueStatus.PENDING,   // job is queued 
            'H': QueueStatus.HOLD,      // job is held
            'S': QueueStatus.HOLD,      // job is suspended 
            'U': QueueStatus.HOLD,      // job is suspended due to workstation becoming busy
            'W': QueueStatus.HOLD,      // job is waiting 
            'T': QueueStatus.HOLD,      // job is in transition
            'M': QueueStatus.HOLD,      // job was moved to another server
    ]

    @Override
    protected QueueStatus decode(String status) {
        DECODE_STATUS.get(status)
    }

}
