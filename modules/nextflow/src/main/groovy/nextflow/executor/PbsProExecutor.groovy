/*
 * Copyright 2013-2024, Seqera Labs
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

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.processor.TaskArrayRun
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
@CompileStatic
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

        if( task instanceof TaskArrayRun ) {
            final arraySize = task.getArraySize()
            result << '-J' << "0-${arraySize - 1}".toString()
        }

        // when multiple competing directives are provided, only the first one will take effect
        // therefore clusterOptions is added as first to give priority over other options as expected
        // by the clusterOptions semantics -- see https://github.com/nextflow-io/nextflow/pull/2036
        addClusterOptionsDirective(task.config, result)

        result << '-N' << getJobNameFor(task)

        result << '-o' << (task.isArray() ? '/dev/null' : quote(task.workDir.resolve(TaskRun.CMD_LOG)))
        result << '-j' << 'oe'

        // the requested queue name
        if( task.config.queue ) {
            result << '-q'  << (String)task.config.queue
        }

        def res = []
        if( task.config.hasCpus() || task.config.getMemory() ) {
            res << "ncpus=${task.config.getCpus()}".toString()
        }
        if( task.config.getMemory() ) {
            // https://www.osc.edu/documentation/knowledge_base/out_of_memory_oom_or_excessive_memory_usage
            res << "mem=${task.config.getMemory().getMega()}mb".toString()
        }
        if( res ) {
            if( matchOptions(task.config.getClusterOptionsAsString()) ) {
                log.warn1 'cpus and memory directives are ignored when clusterOptions contains -l option\ntip: clusterOptions = { "-l select=1:ncpus=${task.cpus}:mem=${task.memory.toMega()}mb:..." }'
            }
            else {
                result << '-l' << "select=1:${res.join(':')}".toString()
            }
        }

        // max task duration
        if( task.config.getTime() ) {
            final duration = task.config.getTime()
            result << "-l" << "walltime=${duration.format('HH:mm:ss')}".toString()
        }

        // add account from config
        final account = session.getExecConfigProp(getName(), 'account', null) as String
        if( account ) {
            result << '-P' << account
        }

        return result
    }

    @Override
    protected List<String> queueStatusCommand(Object queue) {
        String cmd = 'qstat -f ' 
        if( queue ) {
            cmd += queue
        } else {
            cmd += '$( qstat -B | grep -E -v \'(^Server|^---)\' | awk -v ORS=\' \' \'{print \"@\"\$1}\' )'
        }
        return ['bash','-c', "set -o pipefail; $cmd | { grep -E '(Job Id:|job_state =)' || true; }".toString()]
    }

    // see https://www.pbsworks.com/pdfs/PBSRefGuide18.2.pdf
    // table 8.1
    static private Map<String,QueueStatus> DECODE_STATUS = [
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

    @Override
    String getArrayIndexName() {
        return 'PBS_ARRAY_INDEX'
    }

}
