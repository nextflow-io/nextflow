/*
 * Copyright 2013-2024, Seqera Labs
 * Copyright 2022, CEA-CNRGH
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
import nextflow.processor.TaskRun
/**
 * Processor for BRIDGE resource manager (DRAFT)
 *
 * See https://github.com/cea-hpc/bridge 
 *
 *
 * @author Eric Bonnet <eric.d.bonnet@gmail.com>
 */
@Slf4j
@CompileStatic
class BridgeExecutor extends AbstractGridExecutor {

    //  submission pattern example: Submitted Batch Session 1277017
    static private Pattern SUBMIT_REGEX = ~/Submitted Batch Session (\d+)/

    /**
     * Gets the directives to submit the specified task to the cluster for execution
     *
     * @param task A {@link TaskRun} to be submitted
     * @param result The {@link List} instance to which add the job directives
     * @return A {@link List} containing all directive tokens and values.
     */
    protected List<String> getDirectives(TaskRun task, List<String> result) {

        String job_name = ""
        // parenthesis are not compatible with bridge submission commands
        job_name = getJobNameFor(task)
        job_name = job_name.replace("(", "")
        job_name = job_name.replace(")", "")

        result << '-r' << job_name 
        result << '-o' << quote(task.workDir.resolve(TaskRun.CMD_LOG)) 

        // number of cores per parallel task to allocate 
        if( task.config.getCpus() > 1 ) {
            result << '-c' << task.config.getCpus().toString()
        }

        // maximum walltime of the batch job in seconds
        if( task.config.getTime() ) {
            result << '-T' << task.config.getTime().toSeconds().toString() 
        }

        // maximum memory amount required per allocated core in Mo (default is chosen by the underlying system)
        if( task.config.getMemory() ) {
            result << '-M' << task.config.getMemory().toMega().toString() 
        }

        // the requested partition (a.k.a queue) name
        if( task.config.queue ) {
            result << '-q' << (task.config.queue.toString())
        }

        // other cluster options 
        addClusterOptionsDirective(task.config, result)

        return result
    }

    String getHeaderToken() { '#MSUB' }

    /**
     * The command line to submit this job
     *
     * @param task The {@link TaskRun} instance to submit for execution to the cluster
     * @param scriptFile The file containing the job launcher script
     * @return A list representing the submit command line
     */
    @Override
    List<String> getSubmitCommandLine(TaskRun task, Path scriptFile ) {

        ['ccc_msub', scriptFile.getName()]

    }

    /**
     * Parse the string returned by the {@code ccc_msub} command and extract the job ID string
     *
     * @param text The string returned when submitting the job
     * @return The actual job ID string
     */
    @Override
    def parseJobId(String text) {

        for( String line : text.readLines() ) {
            def m = SUBMIT_REGEX.matcher(line)
            if( m.find() ) {
                return m.group(1).toString()
            }
        }

        def id = text.trim()
        if( id.isLong() )
            return id

        throw new IllegalStateException("Invalid BRIDGE submit response:\n$text\n\n")
    }

    @Override
    protected List<String> getKillCommand() { ['ccc_mdel'] }

    @Override
    protected List<String> queueStatusCommand(Object queue) {

        final result = ['ccc_bsstat','-o','batchid,state']

        if( queue )
            result << '-q' << queue.toString()

        final user = System.getProperty('user.name')
        if( user )
            result << '-u' << user
        else
            log.debug "Cannot retrieve current user"

        return result
    }


    /*
     *  Maps job status to nextflow status
     */
    static private Map<String,QueueStatus> STATUS_MAP = [
            'pending': QueueStatus.PENDING, 
            'running': QueueStatus.RUNNING,
            'done': QueueStatus.DONE, 
            'failed': QueueStatus.ERROR, 
            'unknown': QueueStatus.ERROR,
            'suspended': QueueStatus.HOLD,
    ]


    @Override
    protected Map<String, QueueStatus> parseQueueStatus(String text) {

        final result = new LinkedHashMap<String, QueueStatus>()
        if( !text) 
            return result

        for( String line : text.readLines() ) { 

            def cols = line.split(/\s+/)
            if( cols.size() == 2 ) {
                result.put( cols[0], STATUS_MAP.get(cols[1]) )
            }
            else {
                log.debug "[BRIDGE] invalid status line: `$line`"
            }
        }

        return result
    }
}
