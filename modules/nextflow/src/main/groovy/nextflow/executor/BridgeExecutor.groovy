/*
 * Copyright 2020, Seqera Labs
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

        //result << '-D' << quote(task.workDir)
        String job_name = ""
        job_name = getJobNameFor(task)
        job_name = job_name.replace("(", "")
        job_name = job_name.replace(")", "")

        result << '-r' << job_name 
        result << '-o' << quote(task.workDir.resolve(TaskRun.CMD_LOG)) 

        if( task.config.cpus > 1 ) {
            result << '-c' << task.config.cpus.toString()
        }

        if( task.config.time ) {
            //result << '-T' << task.config.getTime().format('HH:mm:ss')
            result << '-T' << task.config.getTime().toSeconds() 
        }

        if( task.config.getMemory() ) {
            //result << '--mem' << task.config.getMemory().toMega().toString() + 'M'
            result << '-M' << task.config.getMemory().toMega().toString() 
        }

        // the requested partition (a.k.a queue) name
        if( task.config.queue ) {
            result << '-q' << (task.config.queue.toString())
        }

        // -- at the end append the command script wrapped file name
        if( task.config.clusterOptions ) {
            result << task.config.clusterOptions.toString() << ''
        }

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
     * Parse the string returned by the {@code sbatch} command and extract the job ID string
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

        final result = ['squeue','--noheader','-o','%i %t', '-t', 'all']

        if( queue )
            result << '-p' << queue.toString()

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
    static private Map STATUS_MAP = [
            'PD': QueueStatus.PENDING,  // (pending)
            'R': QueueStatus.RUNNING,   // (running)
            'CA': QueueStatus.ERROR,    // (cancelled)
            'CF': QueueStatus.PENDING,  // (configuring)
            'CG': QueueStatus.RUNNING,  // (completing)
            'CD': QueueStatus.DONE,     // (completed)
            'F': QueueStatus.ERROR,     // (failed),
            'TO': QueueStatus.ERROR,    // (timeout),
            'NF': QueueStatus.ERROR,    // (node failure)
            'S': QueueStatus.HOLD,      // (job suspended)
            'ST': QueueStatus.HOLD,     // (stopped)
            'PR': QueueStatus.ERROR,    // (Job terminated due to preemption)
            'BF': QueueStatus.ERROR,    // (boot fail, Job terminated due to launch failure)
    ]

    @Override
    protected Map<String, QueueStatus> parseQueueStatus(String text) {

        def result = [:]

        text.eachLine { String line ->
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
