/*
 * Copyright 2022-2025, Seqera Labs
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
 * Processor for Flux Framework executor
 *
 * See https://flux-framework.org
 *
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * @author Vanessa Sochat <sochat1@llnl.gov>
 */
@Slf4j
@CompileStatic
class FluxExecutor extends AbstractGridExecutor {

    static final private Pattern SUBMIT_REGEX = ~/((ƒ|f).+)/

    /**
     * Gets the directives to submit the specified task to the cluster for execution
     *
     * @param task A {@link TaskRun} to be submitted
     * @param result The {@link List} instance to which add the job directives
     * @return A {@link List} containing all directive tokens and values.
     */
    protected List<String> getDirectives(TaskRun task, List<String> result) {
        return result
    }

    // Flux does not require a special token or header
    String getHeaderToken() { null }

    /**
     * The command line to submit this job
     *
     * @param task The {@link TaskRun} instance to submit for execution to the cluster
     * @param scriptFile The file containing the job launcher script
     * @return A list representing the submit command line
     */
    @Override
    List<String> getSubmitCommandLine(TaskRun task, Path scriptFile ) {

        List<String> result = ['flux', 'submit']
        result << '--setattr=cwd=' + quote(task.workDir)
        result << '--job-name="' + getJobNameFor(task) + '"'

        // Only write output to file if user doesn't want written entirely to terminal
        Boolean terminalOutput = session.config.navigate('flux.terminalOutput') as Boolean
        if ( !terminalOutput ) {
            result << '--output=' + quote(task.workDir.resolve(TaskRun.CMD_LOG))  // -o OUTFILE
        }

        if( task.config.getCpus() > 1 ) {
            result << '--cores-per-task=' + task.config.getCpus().toString()
        }

        // Time limit in minutes when no units provided
        if( task.config.getTime() ) {
            result << '--time-limit=' + task.config.getTime().format('mm')
        }

        // Flux does not support asking for specific memory
        if( task.config.getMemory() ) {
            log.debug "Custom memory request is not currently supported by Flux."
        }

        // the requested partition (a.k.a queue) name
        if( task.config.queue ) {
            result << '--queue=' + (task.config.queue.toString())
        }

        // Any extra cluster options the user wants!
        addClusterOptionsDirective(task.config, result)

        result << '/bin/bash' << scriptFile.getName()
        return result
    }

    @Override
    protected void addClusterOptionsDirective(TaskConfig config, List<String> result) {
        final opts = config.getClusterOptions()
        final str = opts instanceof Collection ? opts.join(' ') : opts?.toString()

        if( str ) {
            // Split by space
            for (String item : str.tokenize(' ')) {
                if ( item ) {
                    result << item.stripIndent(true).trim()
                }
            }
        }
    }

    /**
     * Parse the string returned by the {@code flux submit} command and extract the job ID string
     *
     * @param text The string returned when submitting the job
     * @return The actual job ID string
     */
    @Override
    def parseJobId(String text) {

        // Parse the special "f" first
        for( String line : text.readLines() ) {
            def m = SUBMIT_REGEX.matcher(line)
            if( m.find() ) {
                return m.group(1).toString()
            }
        }

        // Fall back to just a jobid
        def id = text.trim()
        if( id.isLong() )
            return id

        throw new IllegalStateException("Invalid Flux submit response:\n$text\n\n")
    }

    @Override
    protected List<String> getKillCommand() { ['flux', 'job', 'cancel'] }

    @Override
    protected List<String> queueStatusCommand(Object queue) {

        // Look at jobs from last 15 minutes
        String command = 'flux jobs --suppress-header --format="{id.f58} {status_abbrev}" --since="-15m"'

        if( queue )
            command += ' --queue=' + queue.toString()

        final user = System.getProperty('user.name')
        if( user )
            command += ' --user=' + user
        else
            log.debug "Cannot retrieve current user"

        final result = ['sh', '-c', command]
        return result
    }

    /*
     *  Maps Flux job status to nextflow status
     *  see https://flux-framework.readthedocs.io/projects/flux-core/en/latest/man1/flux-jobs.html#job-status
     */
    static private Map<String,QueueStatus> STATUS_MAP = [
            'D': QueueStatus.HOLD,      // (depend)
            'R': QueueStatus.RUNNING,   // (running)
            'S': QueueStatus.PENDING,   // (scheduled)
            'C': QueueStatus.DONE,      // (cleanup)
            'I': QueueStatus.DONE,      // (inactive)
    ]

    @Override
    protected Map<String, QueueStatus> parseQueueStatus(String text) {

        final result = new LinkedHashMap<String, QueueStatus>()

        text.eachLine { String line ->
            def cols = line.split(/\s+/)
            if( cols.size() >= 2 ) {
                result.put( cols[0], STATUS_MAP.get(cols[1]) )
            }
            else {
                log.debug "[Flux] invalid status line: `$line`"
            }
        }

        return result
    }
}
