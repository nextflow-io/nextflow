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

import groovy.util.logging.Slf4j
import nextflow.processor.TaskRun
import nextflow.util.Escape

/**
 * Implements a executor for Moab batch scheduler cluster
 *
 * http://www.adaptivecomputing.com/moab-hpc-basic-edition/
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class MoabExecutor extends AbstractGridExecutor {

    /**
     * Gets the directives to submit the specified task to the cluster for execution
     *
     * See http://docs.adaptivecomputing.com/maui/commands/msub.php
     *
     * @param task A {@link TaskRun} to be submitted
     * @param result The {@link List} instance to which add the job directives
     * @return A {@link List} containing all directive tokens and values.
     */
    protected List<String> getDirectives( TaskRun task, List<String> result ) {
        assert result !=null

        result << '-N' << getJobNameFor(task)
        result << '-o' << quote(task.workDir.resolve(TaskRun.CMD_LOG))
        result << '-j' << 'oe'

        // the requested queue name
        if( task.config.queue ) {
            result << '-q'  << (String)task.config.queue
        }

        if( task.config.cpus > 1 ) {
            result << '-l' << "nodes=1:ppn=${task.config.cpus}"
        }

        // max task duration
        if( task.config.time ) {
            final duration = task.config.getTime()
            result << "-l" << "walltime=${duration.format('HH:mm:ss')}"
        }

        // task max memory
        if( task.config.memory ) {
            // https://www.osc.edu/documentation/knowledge_base/out_of_memory_oom_or_excessive_memory_usage
            result << "-l" << "mem=${task.config.memory.toString().replaceAll(/[\s]/,'').toLowerCase()}"
        }

        // -- at the end append the command script wrapped file name
        if( task.config.clusterOptions ) {
            result << task.config.clusterOptions.toString() << ''
        }

        return result
    }

    @Override
    String getHeaders( TaskRun task ) {
        String result = super.getHeaders(task)
        result += "NXF_CHDIR=${Escape.path(task.workDir)}\n"
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
        def cmd = new ArrayList(5)
        cmd << 'msub'
        cmd << '--xml'
        cmd << scriptFile.name
        return cmd
    }

    protected String getHeaderToken() { '#MSUB' }

    /**
     * Parse the string returned by the {@code qsub} command and extract the job ID string
     *
     * @param text The string returned when submitting the job
     * @return The actual job ID string
     */
    @Override
    def parseJobId( String text ) {
        String result
        try {
            result = new XmlSlurper()
                            .parseText(text)
                            .job[0]
                            .@JobID
                            .toString()
        }
        catch (Exception e) {
            throw new IllegalArgumentException("Invalid Moab submit response:\n$text\n\n", e)
        }

        if( !result )
            throw new IllegalArgumentException("Missing Moab submit job ID:\n$text\n\n")

        return result
    }

    @Override
    protected List<String> getKillCommand() { ['mjobctl', '-c'] }

    @Override
    protected List<String> queueStatusCommand(Object queue) {
        return ['showq', '--xml', '-w', "user="+System.getProperty('user.name')]
    }

    static private Map<String,QueueStatus> DECODE_STATUS = [
            'Running': QueueStatus.RUNNING,
            'Starting': QueueStatus.RUNNING,
            'Idle': QueueStatus.HOLD,
            'UserHold': QueueStatus.HOLD,
            'SystemHold': QueueStatus.HOLD,
            'Deferred': QueueStatus.HOLD,
            'NotQueued': QueueStatus.ERROR,
            'BatchHold': QueueStatus.ERROR,
    ]

    protected QueueStatus decode(String status) {
        DECODE_STATUS.get(status)
    }

    @Override
    protected Map<String, QueueStatus> parseQueueStatus(String xmlStatus) {
        // parse XML string and and decode job states
        final result = new LinkedHashMap()
        try {
            def data = new XmlSlurper().parseText(xmlStatus)
            for( def queue : data.queue ) {
                if( !queue.@count?.toInteger() )
                    continue

                def state = queue.@option?.toString()
                if( state=="active")
                    parseQueueJobNodes(queue, result, QueueStatus.RUNNING )
                else if( state=="eligible" )
                    parseQueueJobNodes(queue, result, QueueStatus.PENDING )
                else
                    parseQueueJobNodes(queue, result )
            }
        }
        catch (Exception e) {
            throw e
        }
        return result
    }

    protected void parseQueueJobNodes(queueNode, Map result, QueueStatus state=null) {
        for( def entry : queueNode.job ) {
            final value = state ?: decode(entry.@State?.toString())
            result.put( entry.@JobID?.toString(), value )
        }
    }


    protected List<String> killTaskCommand(def jobId) {
        final result = getKillCommand()
        if( jobId instanceof Collection ) {
            result.add( jobId.join(',') )
            log.trace "Kill command: ${result}"
        }
        else {
            result.add(jobId.toString())
        }
        return result
    }

}
