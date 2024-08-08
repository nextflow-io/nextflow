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

import groovy.transform.CompileStatic
import nextflow.fusion.FusionHelper
import nextflow.processor.TaskArrayRun
import nextflow.processor.TaskRun
/**
 * Execute a task script by running it on the SGE/OGE cluster
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class SgeExecutor extends AbstractGridExecutor implements TaskArrayExecutor {

    /**
     * Gets the directives to submit the specified task to the cluster for execution
     *
     * @param task A {@link TaskRun} to be submitted
     * @param result The {@link List} instance to which add the job directives
     * @return A {@link List} containing all directive tokens and values.
     */
    protected List<String> getDirectives(TaskRun task, List<String> result) {

        if( task instanceof TaskArrayRun ) {
            final arraySize = task.getArraySize()
            result << '-t' << "1-${arraySize}".toString()
        }

        result << '-N' << getJobNameFor(task)

        result << '-o' << (task.isArray() ? '/dev/null' : quote(task.workDir.resolve(TaskRun.CMD_LOG)))
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
            result << "-pe" << "${task.config.penv} ${task.config.getCpus()}".toString()
        }
        else if( task.config.getCpus()>1 ) {
            result << "-l" << "slots=${task.config.getCpus()}".toString()
        }

        // max task duration
        if( task.config.getTime() ) {
            final time = task.config.getTime()
            result << "-l" << "h_rt=${time.format('HH:mm:ss')}".toString()
        }

        // task max memory
        if( task.config.getMemory() ) {
            final mem = "${task.config.getMemory().mega}M".toString()
            result << "-l" << "h_rss=$mem,mem_free=$mem".toString()
        }

        // -- at the end append the command script wrapped file name
        addClusterOptionsDirective(task.config, result)

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
        return pipeLauncherScript()
                ? List.of('qsub', '-')
                : List.of('qsub', '-terse', scriptFile.name)
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

            if( (id=entry.tokenize('.').get(0)).isLong() )
                return id

            if( entry.startsWith('Your job') && entry.endsWith('has been submitted') && (id=entry.tokenize().get(2)) )
                return id

            if( entry.startsWith('Your job array') && entry.endsWith('has been submitted') && (id=entry.tokenize().get(3)) )
                return id.tokenize('.').get(0)
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

    static protected Map<String,QueueStatus> DECODE_STATUS = [
        't': QueueStatus.RUNNING,
        'r': QueueStatus.RUNNING,
        'R': QueueStatus.RUNNING,
        'hr': QueueStatus.RUNNING,
        'qw': QueueStatus.PENDING,
        'h': QueueStatus.PENDING,
        'w': QueueStatus.PENDING,
        'P': QueueStatus.PENDING,
        'N': QueueStatus.PENDING,
        'S': QueueStatus.HOLD,
        's': QueueStatus.HOLD,
        'T': QueueStatus.HOLD,
        'Tr': QueueStatus.HOLD,
        'hqw': QueueStatus.HOLD,
        'Eqw': QueueStatus.ERROR,
        'E': QueueStatus.ERROR
    ]

    @Override
    protected Map<String, QueueStatus> parseQueueStatus(String text) {

        final result = new LinkedHashMap<String, QueueStatus>()
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

    @Override
    protected boolean pipeLauncherScript() {
        return isFusionEnabled()
    }

    @Override
    boolean isFusionEnabled() {
        return FusionHelper.isFusionEnabled(session)
    }

    @Override
    String getArrayIndexName() {
        return 'SGE_TASK_ID'
    }

    @Override
    int getArrayIndexStart() {
        return 1
    }

    @Override
    String getArrayTaskId(String jobId, int index) {
        return "${jobId}.${index}"
    }

}
