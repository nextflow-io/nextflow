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
/*
  Implemented by Satoshi Ohshima <ohshima@cc.kyushu-u.ac.jp>
  I referred to the other executors groovy implementations.
  Resource requests and other job characteristics can be controlled via the following process directives:
  - time
  - clusterOptions
  queue (resource group), cpu, node, and other options should be indicated by clusterOptions.
  This is because they depend on target systems (required options are not the same) and can be controlled by "-L" options in arguments of pjsub command.

  example of nextflow.config on Genkai (Kyushu University)
    process {
      executor = 'tcs'
      time = '00:30:00'
      clusterOptions = '-L rscgrp=a-batch -L vnode-core=4'
    }

  example of nextflow.config on Flow (Nagoya University)
    process {
      executor = 'tcs'
      time = '00:30:00'
      clusterOptions = '-L rscunit=cx -L rscgrp=cx-share -L gpu=1'
    }

 */

package nextflow.executor

import java.nio.file.Path
import java.util.regex.Pattern

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.processor.TaskArrayRun
import nextflow.processor.TaskRun
/*
 * Implements a executor for Fujitsu TCS
 *
 * See https://software.fujitsu.com/jp/manual/manualindex/p21000154.html (in Japanese)
 */
@Slf4j
@CompileStatic
class TcsExecutor extends AbstractGridExecutor implements TaskArrayExecutor {

	static private Pattern SUBMIT_REGEX = ~/\[INFO\] PJM 0000 pjsub Job (\d+) submitted./
    /* modify jobname for TCS on Fugaku */
    static String modName(String name1){
		String name2 = name1.replaceAll("\\(", "")
		String name3 = name2.replaceAll("\\)", "")
		return name3
	}

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
            result << '--bulk --sparam' << "0-${arraySize - 1}".toString()
        }

        result << '-N' << modName(getJobNameFor(task))

		result << '-o' << (task.isArray() ? '/dev/null' : quote(task.workDir.resolve(TaskRun.CMD_LOG)))
        result << '-j' << ''
        result << '-S' << ''

        // max task duration
        if( task.config.getTime() ) {
            final duration = task.config.getTime()
            result << "-L" << "elapse=${duration.format('HH:mm:ss')}".toString()
        }

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
        [ 'pjsub', '-N', modName(getJobNameFor(task)), scriptFile.getName() ]
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
        //String cmd = 'pjstat | grep -v JOB_ID'
        //if( queue ) cmd += ' ' + queue
        //return ['bash','-c', "set -o pipefail; $cmd | { grep -E '(Job Id:|job_state =)' || true; }".toString()]
		final result = ['pjstat']
		return result
    }

    static private Map<String,QueueStatus> DECODE_STATUS = [
            'ACC': QueueStatus.PENDING,
            'QUE': QueueStatus.PENDING,
            'RNA': QueueStatus.PENDING,
            'RUN': QueueStatus.RUNNING,
            'RNO': QueueStatus.RUNNING,
            'EXT': QueueStatus.RUNNING,
            'CCL': QueueStatus.DONE,
            'HLD': QueueStatus.HOLD,
            'ERR': QueueStatus.ERROR
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
        value ? SUBMIT_REGEX.matcher(value).find() : null
    }

    @Override
    String getArrayIndexName() {
        return 'TCS_SUBJOBID'
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
