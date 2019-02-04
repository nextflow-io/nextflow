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

import groovy.transform.InheritConstructors
import nextflow.processor.TaskRun
/**
 * HTCondor executor
 *
 * See https://research.cs.wisc.edu/htcondor/
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CondorExecutor extends AbstractGridExecutor {

    static final public String CMD_CONDOR = '.command.condor'

    final protected BashWrapperBuilder createBashWrapperBuilder(TaskRun task) {
        // creates the wrapper script
        final builder = new CondorWrapperBuilder(task)
        builder.manifest = getDirectivesText(task)
        return builder
    }

    protected String getDirectivesText(TaskRun task) {
        def lines = getDirectives(task)
        lines << ''
        lines.join('\n')
    }

    @Override
    protected String getHeaderToken() {
        throw new UnsupportedOperationException()
    }

    @Override
    protected List<String> getDirectives(TaskRun task, List<String> result) {

        result << "universe = vanilla"
        result << "executable = ${TaskRun.CMD_RUN}"
        result << "log = ${TaskRun.CMD_LOG}"
        result << "getenv = true"

        if( task.config.getCpus()>1 ) {
            result << "request_cpus = ${task.config.getCpus()}"
            result << "machine_count = 1"
        }

        if( task.config.getMemory() ) {
            result << "request_memory = ${task.config.getMemory()}"
        }

        if( task.config.getDisk() ) {
            result << "request_disk = ${task.config.getDisk()}"
        }

        if( task.config.getTime() ) {
            result << "periodic_remove = (RemoteWallClockTime - CumulativeSuspensionTime) > ${task.config.getTime().toSeconds()}"
        }

        if( task.config.clusterOptions ) {
            def opts = task.config.clusterOptions
            if( opts instanceof Collection ) {
                result.addAll(opts)
            }
            else {
                result.addAll( opts.toString().tokenize(';\n').collect{ it.trim() })
            }
        }

        result<< "queue"

    }

    @Override
    List<String> getSubmitCommandLine(TaskRun task, Path scriptFile) {
        return ['condor_submit', '--terse', CMD_CONDOR]
    }

    @Override
    def parseJobId(String text) {
        text.tokenize(' -')[0]
    }

    @Override
    protected List<String> getKillCommand() {
        ['condor_rm']
    }

    @Override
    protected List<String> queueStatusCommand(Object queue) {
        ["condor_q"]
    }


    static protected Map DECODE_STATUS = [
            'U': QueueStatus.PENDING,   // Unexpanded
            'I': QueueStatus.PENDING,   // Idle
            'R': QueueStatus.RUNNING,   // Running
            'X': QueueStatus.ERROR,     // Removed
            'C': QueueStatus.DONE,      // Completed
            'H': QueueStatus.HOLD,      // Held
            'E': QueueStatus.ERROR      // Error
    ]


    @Override
    protected Map<String, QueueStatus> parseQueueStatus(String text) {
        def result = [:]
        if( !text ) return result

        boolean started = false
        def itr = text.readLines().iterator()
        while( itr.hasNext() ) {
            String line = itr.next()
            if( !started ) {
                started = line.startsWith(' ID ')
                continue
            }

            if( !line.trim() ) {
                break
            }

            def cols = line.tokenize(' ')
            def id = cols[0]
            def st = cols[5]
            result[id] = DECODE_STATUS[st]
        }

        return result
    }


    @InheritConstructors
    static class CondorWrapperBuilder extends BashWrapperBuilder {

        String manifest

        Path build() {
            final wrapper = super.build()
            // give execute permission to wrapper file
            wrapper.setExecutable(true)
            // save the condor manifest
            workDir.resolve(CMD_CONDOR).text = manifest
            return wrapper
        }

    }
}
