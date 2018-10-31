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

        result << '-N' << getJobNameFor(task)
        result << '-o' << quote(task.workDir.resolve(TaskRun.CMD_LOG))
        result << '-j' << 'oe'

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

    @Override
    String getHeaders( TaskRun task ) {
        String result = super.getHeaders(task)
        result += "cd ${quote(task.workDir)}\n"
        return result
    }

    @Override
    String getJobNameFor( TaskRun task ) {
        def result = super.getJobNameFor(task)
        // some implementations do not allow parenthesis in the job name -- see #271
        result = result.replace('(','').replace(')','')
        // PBS does not allow more than 15 characters for the job name string
        result && result.size()>15 ? result.substring(0,15) : result
    }
    /**
     * The command line to submit this job
     *
     * @param task The {@link TaskRun} instance to submit for execution to the cluster
     * @param scriptFile The file containing the job launcher script
     * @return A list representing the submit command line
     */
    List<String> getSubmitCommandLine(TaskRun task, Path scriptFile ) {
        // in some PBS implementation the submit command will fail if the script name starts with a dot eg `.command.run`
        // add the `-N <job name>` to fix this -- see issue #228
        [ 'qsub', '-N', getJobNameFor(task), scriptFile.getName() ]
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
        if( result && !result.contains('\n')) {
            return result
        }

        throw new IllegalArgumentException("Invalid PBS/Torque submit response:\n$text\n\n")
    }

    @Override
    protected List<String> getKillCommand() { ['qdel'] }

    @Override
    protected List<String> queueStatusCommand(Object queue) {
        String cmd = 'qstat -f -1'
        if( queue ) cmd += ' ' + queue
        return ['bash','-c', "set -o pipefail; $cmd | { egrep '(Job Id:|job_state =)' || true; }".toString()]
    }

    static private Map DECODE_STATUS = [
            'C': QueueStatus.DONE,
            'R': QueueStatus.RUNNING,
            'Q': QueueStatus.PENDING,
            'H': QueueStatus.HOLD,
            'S': QueueStatus.HOLD
    ]

    @Override
    protected Map<String, QueueStatus> parseQueueStatus(String text) {

        final JOB_ID = 'Job Id:'
        final JOB_STATUS = 'job_state ='
        final result = [:]

        String id = null
        String status = null
        text.eachLine { line ->
            if( line.startsWith(JOB_ID) ) {
                id = fetchValue(JOB_ID, line)
            }
            else if( id ) {
                status = fetchValue(JOB_STATUS, line)
            }
            result.put( id, DECODE_STATUS[status] ?: AbstractGridExecutor.QueueStatus.UNKNOWN )
        }

        return result
    }

    static String fetchValue( String prefix, String line ) {
        final p = line.indexOf(prefix)
        return p!=-1 ? line.substring(p+prefix.size()).trim() : null
    }

}
