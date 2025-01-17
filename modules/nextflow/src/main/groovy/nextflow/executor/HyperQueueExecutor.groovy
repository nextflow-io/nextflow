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
 *
 */

package nextflow.executor

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.processor.TaskRun
import nextflow.util.ServiceName
import nextflow.util.TupleHelper
/**
 * Implements executor for HyperQueue
 *  https://github.com/It4innovations/hyperqueue/
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * @author Henrik Nortamo <henrik.nortamo@csc.fi>
 */
@CompileStatic
@ServiceName('hq')
@Slf4j
class HyperQueueExecutor extends AbstractGridExecutor {

    private static final STATUSES = [
        'CANCELED': QueueStatus.DONE,
        'FAILED': QueueStatus.ERROR,
        'FINISHED': QueueStatus.DONE,
        'RUNNING': QueueStatus.RUNNING,
        'WAITING': QueueStatus.HOLD
    ]

    @Override
    protected String getHeaderToken() {
        return '#HQ'
    }

    @Override
    protected void register() {
        super.register()
        log.warn 'The support for HyperQueue is an experimental feature and it may change in a future release'
    }

    @Override
    protected List<String> getDirectives(TaskRun task, List<String> result) {

        result << '--name' << getJobNameFor(task)
        result << '--log' << quote(task.workDir.resolve(TaskRun.CMD_LOG))
        result << '--cwd' << quote(task.workDir)

        // No enforcement, Hq just makes sure that the allocated value is below the limit
        if( task.config.getMemory() )
            result << '--resource' << "mem=${task.config.getMemory().toMega()}".toString()
        if( task.config.hasCpus() )
            result << '--cpus'<< task.config.getCpus().toString()
        if( task.config.getTime() )
            result << '--time-limit' << (task.config.getTime().toSeconds() + 'sec')
        if( task.config.getAccelerator() )
            result << '--resource' << "gpus=${task.config.getAccelerator().limit}".toString()

        // -- At the end append the command script wrapped file name
        addClusterOptionsDirective(task.config, result)

        return result
    }

    @Override
    List<String> getSubmitCommandLine(TaskRun task, Path scriptFile) {
        return TupleHelper.listOf('hq', '--output-mode=quiet', 'submit', '--directives=file', scriptFile.getName())
    }

    @Override
    def parseJobId(String text) {
        try {
            text as Integer as String
        }
        catch (Exception e) {
            throw new IllegalArgumentException("HyperQueue submit failed or response is invalid:\n$text\n\n")
        }
    }

    @Override
    protected List<String> getKillCommand() {
        return TupleHelper.listOf('hq', 'job', 'cancel')
    }

    @Override
    protected List<String> killTaskCommand(def jobId) {
        final result = getKillCommand()
        if( jobId instanceof Collection ) {
            result.add(jobId.join(','))
        }
        else {
            result.add(jobId.toString())
        }
        return result
    }

    @Override
    protected List<String> queueStatusCommand(Object queue) {
        return TupleHelper.listOf('hq', '--output-mode=quiet', 'job', 'list', '--all')
    }

    @Override
    protected Map<String, QueueStatus> parseQueueStatus(String text) {
        final result = new HashMap()
        for( String line : text.split('\n') ) {
            final tokens = line.split(' ')
            final id = tokens[0]
            final status = STATUSES[tokens[1]]
            result.put(id, status)
        }
        return result
    }

}
