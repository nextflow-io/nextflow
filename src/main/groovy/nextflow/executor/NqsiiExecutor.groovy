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
import java.nio.file.Path

import nextflow.processor.TaskRun
/**
 * Execute a task script by running it on the NQSII cluster
 *
 * Read more 
 *   https://www.rz.uni-kiel.de/en/our-portfolio/hiperf/nec-linux-cluster
 *   https://wickie.hlrs.de/platforms/index.php/Batch_system
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * @author Till Bayer <till.bayer@gmail.com>
 */
class NqsiiExecutor extends AbstractGridExecutor {

    /**
     * Gets the directives to submit the specified task to the cluster for execution
     *
     * @param task A {@link TaskRun} to be submitted
     * @param result The {@link List} instance to which add the job directives
     * @return A {@link List} containing all directive tokens and values.
     */
    protected List<String> getDirectives(TaskRun task, List<String> result) {

        result << '-N' << getJobNameFor(task)
        result << '-o' << quote(task.workDir.resolve(TaskRun.CMD_LOG))
        result << '-j' << 'o'  // stderr to stdin
        result << '-b' << '1'    // set number of nodes to 1

        // the requested queue name
        if( task.config.queue ) {
            result << '-q' << (task.config.queue as String)
        }

        //number of cpus for multiprocessing/multi-threading
        if( task.config.cpus>1 ) {
            result << "-l" << "cpunum_job=${task.config.cpus}"
        } else {
            result << "-l" << "cpunum_job=1"
        }

        // max task duration
        if( task.config.time ) {
            final time = task.config.getTime()
            result << "-l" << "elapstim_req=${time.format('HH:mm:ss')}"
        }

        // task max memory
        if( task.config.memory ) {
            result << "-l" << "memsz_job=${task.config.memory.toString().replaceAll(/[\s]/,'').toLowerCase()}"
        }

        // -- at the end append the command script wrapped file name
        if( task.config.clusterOptions ) {
            result << task.config.clusterOptions.toString() << ''
        }

        return result
    }


    /*
     * Prepare the 'qsub' cmdline
     */
    List<String> getSubmitCommandLine(TaskRun task, Path scriptFile ) {
        return ['qsub', scriptFile.name]
    }

    protected String getHeaderToken() { '#PBS' }

    @Override
    String getHeaders( TaskRun task ) {
        String result = super.getHeaders(task)
        result += "cd ${quote(task.workDir)}\n"
        return result
    }

    @Override
    String getJobNameFor( TaskRun task ) {
        def result = super.getJobNameFor(task)
        // NQSII does not allow more than 63 characters for the job name string
        result && result.size()>63 ? result.substring(0,63) : result
    }

    /**
     * Parse the string returned by the {@code qsub} command and extract the job ID string
     *
     * @param text The string returned when submitting the job
     * @return The actual job ID string
     */
    @Override
    def parseJobId( String text ) {
    def pattern = ~/Request (\d+).+ submitted to queue.+/
    for( String line : text.readLines() ) {
            def m = pattern.matcher(line)
            if( m.find() ) {
                return m.group(1)
            }
        }

        throw new IllegalStateException("Invalid NQSII submit response:\n$text\n\n")
    }

    @Override
    protected List<String> getKillCommand() { ['qdel'] }

    @Override
    protected List<String> queueStatusCommand(Object queue) {
        def result = ['qstat']
        if( queue )
            result << '-q' << queue.toString()

        return result
    }

    static protected Map DECODE_STATUS = [
            'PRR': QueueStatus.RUNNING,
            'RUN': QueueStatus.RUNNING,
            'STG': QueueStatus.RUNNING,
            'CHK': QueueStatus.RUNNING,
            'EXT': QueueStatus.RUNNING,
            'POR': QueueStatus.RUNNING,
            'ARI': QueueStatus.PENDING,
            'FWD': QueueStatus.PENDING,
            'MIG': QueueStatus.PENDING,
            'QUE': QueueStatus.PENDING,
            'WAT': QueueStatus.PENDING,
            'GQD': QueueStatus.PENDING,
            'SUS': QueueStatus.HOLD,
            'HLD': QueueStatus.HOLD,
            'HOL': QueueStatus.HOLD
    ]

    @Override
    protected Map<String, QueueStatus> parseQueueStatus(String text) {

        def result = [:]
        text?.eachLine{ String row, int index ->
            if( index< 2 ) return
            def cols = row.trim().split(/\s+/)
            if( cols.size()>6 ) {
            def id = cols[0].tokenize('.')
                result.put( id[0], DECODE_STATUS[cols[5]] )
            }
        }

        return result
    }

}
