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

import java.nio.file.Path
import java.util.regex.Pattern

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.processor.TaskArrayRun
import nextflow.processor.TaskRun
/**
 * Implements a executor for PBS/Torque cluster
 *
 * See http://www.pbsworks.com
 */
@Slf4j
@CompileStatic
class PbsExecutor extends AbstractGridExecutor implements TaskArrayExecutor {

    private static Pattern OPTS_REGEX = ~/(?:^|\s)-l.+/

    /**
     * Gets the directives to submit the specified task to the cluster for execution
     *
     * @param task A {@link TaskRun} to be submitted
     * @param result The {@link List} instance to which add the job directives
     * @return A {@link List} containing all directive tokens and values.
     */
    protected List<String> getDirectives( TaskRun task, List<String> result ) {
        assert result !=null

        if( task instanceof TaskArrayRun ) {
            final arraySize = task.getArraySize()
            result << '-t' << "0-${arraySize - 1}".toString()
        }

        result << '-N' << getJobNameFor(task)

        result << '-o' << (task.isArray() ? '/dev/null' : quote(task.workDir.resolve(TaskRun.CMD_LOG)))
        result << '-j' << 'oe'

        // the requested queue name
        if( task.config.queue ) {
            result << '-q'  << (String)task.config.queue
        }

        // task cpus
        if( task.config.getCpus() > 1 ) {
            if( matchOptions(task.config.getClusterOptionsAsString()) ) {
                log.warn1 'cpus directive is ignored when clusterOptions contains -l option\ntip: clusterOptions = { "-l nodes=1:ppn=${task.cpus}:..." }'
            }
            else {
                result << '-l' << "nodes=1:ppn=${task.config.getCpus()}".toString()
            }
        }

        // max task duration
        if( task.config.getTime() ) {
            final duration = task.config.getTime()
            result << "-l" << "walltime=${duration.format('HH:mm:ss')}".toString()
        }

        // task max memory
        if( task.config.getMemory() ) {
            // https://www.osc.edu/documentation/knowledge_base/out_of_memory_oom_or_excessive_memory_usage
            result << "-l" << "mem=${task.config.getMemory().toString().replaceAll(/[\s]/,'').toLowerCase()}".toString()
        }

        // add account from config
        final account = session.getExecConfigProp(getName(), 'account', null) as String
        if( account ) {
            result << '-P' << account
        }

        // -- at the end append the command script wrapped file name
        addClusterOptionsDirective(task.config, result)

        return result
    }

    @Override
    String sanitizeJobName( String name ) {
        // some implementations do not allow parenthesis in the job name -- see #271
        name = name.replace('(','').replace(')','')
        // PBS does not allow more than 15 characters for the job name string
        name.size()>15 ? name.substring(0,15) : name
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
        return ['bash','-c', "set -o pipefail; $cmd | { grep -E '(Job Id:|job_state =)' || true; }".toString()]
    }

    static private Map<String,QueueStatus> DECODE_STATUS = [
            'C': QueueStatus.DONE,
            'R': QueueStatus.RUNNING,
            'Q': QueueStatus.PENDING,
            'H': QueueStatus.HOLD,
            'S': QueueStatus.HOLD
    ]

    protected QueueStatus decode(String status) {
        DECODE_STATUS.get(status)
    }

    @Override
    protected Map<String, QueueStatus> parseQueueStatus(String text) {

        final JOB_ID = 'Job Id:'
        final JOB_STATUS = 'job_state ='
        final result = new LinkedHashMap<String, QueueStatus>()

        String id = null
        String status = null
        text.eachLine { line ->
            if( line.startsWith(JOB_ID) ) {
                id = fetchValue(JOB_ID, line)
            }
            else if( id ) {
                status = fetchValue(JOB_STATUS, line)
            }
            result.put( id, decode(status) ?: QueueStatus.UNKNOWN )
        }

        return result
    }

    static String fetchValue( String prefix, String line ) {
        final p = line.indexOf(prefix)
        return p!=-1 ? line.substring(p+prefix.size()).trim() : null
    }

    static protected boolean matchOptions(String value) {
        value ? OPTS_REGEX.matcher(value).find() : null
    }

    @Override
    String getArrayIndexName() {
        return 'PBS_ARRAYID'
    }

    @Override
    int getArrayIndexStart() {
        return 0
    }

    @Override
    String getArrayTaskId(String jobId, int index) {
        assert jobId, "Missing 'jobId' argument"
        return jobId.replace('[]', "[$index]")
    }

}
