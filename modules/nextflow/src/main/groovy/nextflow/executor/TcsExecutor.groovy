/*
 * Copyright 2013-2026, Seqera Labs
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

/*
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

/*
 * Implements a executor for Fujitsu Technical Computing Suite
 *
 * https://software.fujitsu.com/jp/manual/manualindex/p21000155e.html
 *
 * @author Satoshi Ohshima <ohshima@cc.kyushu-u.ac.jp>
 */
@Slf4j
@CompileStatic
class TcsExecutor extends AbstractGridExecutor implements TaskArrayExecutor {

    static private final Pattern SUBMIT_REGEX = ~/\[INFO\] PJM 0000 pjsub Job (\d+) submitted./

    /**
     * Modify job name for TCS on Fugaku.
     *
     * @param name
     */
    static String modName(String name) {
        return name
            .replaceAll("\\(", "")
            .replaceAll("\\)", "")
    }

    /**
     * Gets the directives to submit the specified task to the cluster for execution
     *
     * @param task A {@link TaskRun} to be submitted
     * @param result The {@link List} instance to which add the job directives
     * @return A {@link List} containing all directive tokens and values.
     */
    protected List<String> getDirectives(TaskRun task, List<String> result) {
        assert result != null, "Argument 'result' cannot be null"

        result << '-N' << modName(getJobNameFor(task))

        // max task duration
        if( task.config.getTime() ) {
            final duration = task.config.getTime()
            result << "-L" << "elapse=${duration.format('HH:mm:ss')}".toString()
        }

        // output file
        if( task.isArray() ) {
            // If task is array job (bulk job), don't indicate output file (use default setting).
            // * TCS doesn't support /dev/null for output. (depend on system?)
        } else {
            result << '-o' << quote(task.workDir.resolve(TaskRun.CMD_LOG))
        }

        // output options
        result << '-j' << '' // marge stderr to stdout
        result << '-S' << '' // output information

        // -- at the end append the command script wrapped file name
        addClusterOptionsDirective(task.config, result)

        return result
    }

    /**
     * The command line to submit this job
     *
     * @param task The {@link TaskRun} instance to submit for execution to the cluster
     * @param scriptFile The file containing the job launcher script
     * @return A list representing the submit command line
     */
    List<String> getSubmitCommandLine(TaskRun task, Path scriptFile ) {
        if( task instanceof TaskArrayRun ) {
            final arraySize = task.getArraySize()
            return pipeLauncherScript()
                ? List.of('pjsub', '--bulk --sparam ', "0-${arraySize - 1}".toString())
                : List.of('pjsub', scriptFile.getName())
        } else {
            return pipeLauncherScript()
                ? List.of('pjsub')
                : List.of('pjsub', scriptFile.getName())
        }
    }

    protected String getHeaderToken() { '#PJM' }

    /**
     * Parse the string returned by the {@code pjsub} command and extract the job ID string
     *
     * @param text The string returned when submitting the job
     * @return The actual job ID string
     */
    @Override
    def parseJobId( String text ) {
        for( String line : text.readLines() ) {
            log.warn1 line
            def m = SUBMIT_REGEX.matcher(line)
            if( m.find() ) {
                return m.group(1).toString()
            }
        }
        throw new IllegalArgumentException("Invalid TCS submit response:\n$text\n\n")
    }

    @Override
    protected List<String> getKillCommand() { ['pjdel'] }

    @Override
    protected List<String> queueStatusCommand(Object queue) {
        return ['pjstat', '-E']
    }

    static private Map<String,QueueStatus> STATUS_MAP = [
            'ACC': QueueStatus.PENDING, // accepted
            'QUE': QueueStatus.PENDING, // wait for running
            'RNA': QueueStatus.RUNNING, // preparing
            'RUN': QueueStatus.RUNNING, // running
            'RNO': QueueStatus.RUNNING, // cleanup
            'EXT': QueueStatus.DONE,    // finished
            'CCL': QueueStatus.DONE,    // canceled
            'HLD': QueueStatus.HOLD,    // holding
            'ERR': QueueStatus.ERROR,   // error
    ]

    @Override
    protected Map<String, QueueStatus> parseQueueStatus(String text) {
        final result = new LinkedHashMap<String, QueueStatus>()
        text.eachLine { String line ->
            def cols = line.split(/\s+/)
            if( cols.size() > 1 ) {
                result.put( cols[0], STATUS_MAP.get(cols[3]) )
            }
            else {
                log.debug "[TCS] invalid status line: `$line`"
            }
        }

        return result
    }

    @Override
    String getArrayIndexName() {
        return 'PJM_BULKNUM'
    }

    @Override
    int getArrayIndexStart() {
        return 0
    }

    @Override
    String getArrayTaskId(String jobId, int index) {
        assert jobId, "Missing 'jobId' argument"
        return "${jobId}[${index}]"
    }

}
