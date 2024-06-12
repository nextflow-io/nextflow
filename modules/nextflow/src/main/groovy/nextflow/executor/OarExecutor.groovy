/*
 * Copyright 2019, Hospices Civils de Lyon (HCL)
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
import nextflow.processor.TaskConfig
import nextflow.processor.TaskRun
/**
 * Processor for OAR resource manager
 *
 * See https://oar.imag.fr/
 *
 * @author Maxime Vallée <maxime.vallee@chu-lyon.fr>
 *
 */
@Slf4j
@CompileStatic
class OarExecutor extends AbstractGridExecutor {

    /**
     * Gets the directives to submit the specified task to the cluster for execution
     *
     * @param task A {@link TaskRun} to be submitted
     * @param result The {@link List} instance to which add the job directives
     * @return A {@link List} containing all directive tokens and values.
     */
    protected List<String> getDirectives(TaskRun task, List<String> result) {

        result << '-d' << quote(task.workDir)
        result << '-n' << getJobNameFor(task)
        // see discussion https://github.com/nextflow-io/nextflow/issues/1761
        result << '-O' << quote(task.workDir.resolve(TaskRun.CMD_LOG))
        result << '-E' << quote(task.workDir.resolve(TaskRun.CMD_LOG))

        if( task.config.getMemory() ) {
            if( task.config.getMemory().toMega() < 1024 ) {
                result << "-p" << "\"memnode=1\""
            }
            else {
                result << "-p" << "\"memnode=${task.config.getMemory().toGiga().toString()}\"".toString()
            }
        }

        if( task.config.getCpus() > 1) {
            if( task.config.getTime() ) {
                // cpu + time set
                result << "-l" << "/nodes=1/core=${task.config.getCpus().toString()},walltime=${task.config.getTime().format('HH:mm:ss')}".toString()
            }
            else {
                // just cpu set
                result << "-l" << "/nodes=1/core=${task.config.getCpus().toString()}".toString()
            }
        }
        else {
            if( task.config.getTime() ) {
                // just time set
                result << "-l" << "walltime=${task.config.getTime().format('HH:mm:ss')}".toString()
            }
        }
        
        // the requested queue name
        if( task.config.queue ) {
            result << '-q' << (task.config.queue.toString())
        }

        // -- at the end append the command script wrapped file name
        addClusterOptionsDirective(task.config, result)

        return result
    }

    @Override
    void addClusterOptionsDirective(TaskConfig config, List result) {
        final opts = config.getClusterOptions()
        if( opts instanceof Collection ) {
            for( String it in opts ) {
                result << it << ''
            }
        }
        // Options need to be semicolon ";" separated, if several are needed
        else if( opts instanceof CharSequence ) {
            for (String item : opts.toString().tokenize(';')) {
                result << item << ''
            }
        }
    }

    String getHeaderToken() { '#OAR' }

    /**
     * The command line to submit this job
     *
     * @param task The {@link TaskRun} instance to submit for execution to the cluster
     * @param scriptFile The file containing the job launcher script
     * @return A list representing the submit command line
     */
    @Override
    List<String> getSubmitCommandLine(TaskRun task, Path scriptFile ) {
        // Scripts need to be executable
        scriptFile.setPermissions(7,0,0)
        return ["oarsub", "-S", "./${scriptFile.getName()}".toString()]
    }

    /**
     * Parse the string returned by the {@code oarsub} command and extract the job ID string
     *
     * @param text The string returned when submitting the job
     * @return The actual job ID string
     */
    static private Pattern SUBMIT_REGEX = ~/OAR_JOB_ID=(\d+)/
    
    @Override
    def parseJobId(String text) {
        for( String line : text.readLines() ) {
            def m = SUBMIT_REGEX.matcher(line)
            if( m.matches() ) {
                return m.group(1).toString()
            }
        }
        throw new IllegalStateException("Invalid OAR submit response:\n$text\n\n")
    }

    @Override
    protected List<String> getKillCommand() { ['oardel'] }

    @Override
    protected List<String> queueStatusCommand(Object queue) {
        // To have a parsable list of jobs in queue by user
        // see page 21 http://oar.imag.fr/docs/2.5/OAR-Documentation.pdf
        String cmd = 'oarstat -f'
        if( queue ) cmd += ' ' + queue
        return ['sh','-c', "$cmd | grep -E '(Job_Id:|state =)' ".toString()]
    }

    /*
     *  Maps OAR job status to nextflow status
     *  see page 134 http://oar.imag.fr/docs/2.5/OAR-Documentation.pdf
     */
    static private Map<String,QueueStatus> DECODE_STATUS = [
            'toLaunch': QueueStatus.RUNNING,			// (the OAR scheduler has attributed some nodes to the job. So it will be launched.)
            'Launching': QueueStatus.RUNNING,			// (OAR has launched the job and will execute the user command on the first node.)
            'Running': QueueStatus.RUNNING,				// (the user command is executing on the first node.)
            'Finishing': QueueStatus.RUNNING,			// (the user command has terminated and OAR is doing work internally)
            'Waiting': QueueStatus.PENDING,				// (the job is waiting OAR scheduler decision.)
            'toAckReservation': QueueStatus.PENDING,	// (the OAR scheduler must say “YES” or “NO” to the waiting oarsub command because it requested a reservation.)
            'Hold': QueueStatus.HOLD,					// (user or administrator wants to hold the job (oarhold command). So it will not be scheduled by the system.)
            'Suspended': QueueStatus.HOLD,				// (the job was in Running state and there was a request (oarhold with “-r” option) to suspend this job. In this state other jobs can be scheduled on the same resources (these resources has the “suspended_jobs” field to “YES”).)
            'Error': QueueStatus.ERROR,					// (a problem has occurred.)
            'toError': QueueStatus.ERROR,				// (something wrong occurred and the job is going into the error state.)
            'Terminated': QueueStatus.DONE,				// (the job has terminated normally.)
    ]

    protected QueueStatus decode(String status) {
        DECODE_STATUS.get(status)
    }

    @Override
    protected Map<String, QueueStatus> parseQueueStatus(String text) {

        final JOB_ID = 'Job_Id:'
        final JOB_STATUS = 'state ='
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

    @Override
    protected String sanitizeJobName(String name) {
        name.size() > 100 ? name.substring(0,100) : name
    }
}
