/*
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

import nextflow.processor.TaskRun
/**
 * Execute a task script by running it on the SGE/OGE cluster
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class SgeExecutor extends AbstractGridExecutor {

    /**
     * Gets the directives to submit the specified task to the cluster for execution
     *
     * @param task A {@link TaskRun} to be submitted
     * @param result The {@link List} instance to which add the job directives
     * @return A {@link List} containing all directive tokens and values.
     */
    protected List<String> getDirectives(TaskRun task, List<String> result) {

        result << '-wd' << quote(task.workDir)
        result << '-N' << getJobNameFor(task)
        result << '-o' << quote(task.workDir.resolve(TaskRun.CMD_LOG))
        result << '-j' << 'y'
        result << '-terse' << ''    // note: directive need to be returned as pairs

        /*
         * By using command line option -notify SIGUSR1 will be sent to your script prior to SIGSTOP
         * and SIGUSR2 will be sent to your script prior to SIGKILL
         */
        result << '-notify' << ''

        // the requested queue name
        if( task.config.queue ) {
            result << '-q' << (task.config.queue as String)
        }

        //number of cpus for multiprocessing/multi-threading
        if ( task.config.penv ) {
            result << "-pe" << "${task.config.penv} ${task.config.cpus}"
        }
        else if( task.config.cpus>1 ) {
            result << "-l" << "slots=${task.config.cpus}"
        }

        // max task duration
        if( task.config.time ) {
            final time = task.config.getTime()
            result << "-l" << "h_rt=${time.format('HH:mm:ss')}"
        }

        // task max memory
        if( task.config.memory ) {
            final mem = "${task.config.getMemory().mega}M"
            result << "-l" << "h_rss=$mem,mem_free=$mem"
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

        // The '-terse' command line control the output of the qsub command line, when
        // used it only return the ID of the submitted job.
        // NOTE: In some SGE implementations the '-terse' only works on the qsub command line
        // and it is ignored when used in the script job as directive, fir this reason it
        // should not be remove from here
        return ['qsub', '-terse', scriptFile.name]
    }

    protected String getHeaderToken() { '#$' }


    /**
     * Parse the string returned by the {@code qsub} command and extract the job ID string
     *
     * @param text The string returned when submitting the job
     * @return The actual job ID string
     */
    @Override
    def parseJobId( String text ) {
        // return always the last line
        String id
        def lines = text.trim().readLines()
        def entry = lines[-1].trim()
        if( entry ) {
            if( entry.toString().isLong() )
                return entry

            if( entry.startsWith('Your job') && entry.endsWith('has been submitted') && (id=entry.tokenize().get(2)) )
                return id
        }

        throw new IllegalStateException("Invalid SGE submit response:\n$text\n\n")
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
            'r': QueueStatus.RUNNING,
            'qw': QueueStatus.PENDING,
            'hqw': QueueStatus.HOLD,
            'Eqw': QueueStatus.ERROR
    ]

    @Override
    protected Map<String, QueueStatus> parseQueueStatus(String text) {

        def result = [:]
        text?.eachLine{ String row, int index ->
            if( index< 2 ) return
            def cols = row.trim().split(/\s+/)
            if( cols.size()>5 ) {
                result.put( cols[0], DECODE_STATUS[cols[4]] )
            }
        }

        return result
    }

    @Override
    String quote(Path path) {
        // note: SGE does not recognize `\` escape character in the
        // in the path defined as `#$` directives
        // just double-quote paths containing blanks
        def str = path.toString()
        str.indexOf(' ') != -1 ? "\"$str\"" : str
    }


}
