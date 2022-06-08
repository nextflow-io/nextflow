/*
 * Copyright 2020-2022, Seqera Labs
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

import groovy.json.JsonSlurper
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

    private final JsonSlurper parser = new JsonSlurper()

    @Override
    protected String getHeaderToken() {
        return '#HQ'
    }

    @Override
    protected void register() {
        super.register()
        log.warn "The support for HyperQueue is an experimental feature and it may change in a future release"
    }

    @Override
    protected List<String> getDirectives(TaskRun task, List<String> result) {

        result << '--name' << getJobNameFor(task)
        result << '--log' << quote(task.workDir.resolve(TaskRun.CMD_LOG))
        result << '--cwd' << quote(task.workDir)

        // No enforcement, Hq just makes sure that the allocated value is below the limit
        if( task.config.getMemory() )
            result << '--resource' << "mem=${task.config.getMemory().toBytes()}".toString()
        if( task.config.hasCpus() )
            result << '--cpus'<< task.config.getCpus().toString()
        if( task.config.getTime() )
            result << '--time-limit' << (task.config.getTime().toSeconds() + 'sec')
        if( task.config.getAccelerator() )
            result << '--resource' << "gpus=${task.config.getAccelerator().limit}".toString()
        
        // -- At the end append the command script wrapped file name
        if( task.config.clusterOptions ) {
            result << task.config.clusterOptions.toString() << ''
        }

        return result
    }

    @Override
    List<String> getSubmitCommandLine(TaskRun task, Path scriptFile) {
        return TupleHelper.listOf('hq', '--output-mode=json','submit', '--directives=file', scriptFile.getName())
    }

    @Override
    def parseJobId(String text) {
        Map result = null
        try {
            result = parser.parseText(text) as Map
        }
        catch (Exception e) {
            throw new IllegalArgumentException("Invalid HyperQueue submit JSON response:\n$text\n\n")
        }

        if( result.id )
            return result.id as String

        final err = result.error
                    ? "HyperQueue submit error: ${result.error}"
                    : "Invalid HyperQueue submit response:\n$text\n\n"
        throw new IllegalArgumentException(err)
    }

    @Override
    protected List<String> getKillCommand() {
        return TupleHelper.listOf('hq','job', 'cancel')
    }

    // JobIds to the cancel command need to be separated by a comma 
    @Override
    protected List<String> killTaskCommand(def jobId) { 
        final result = getKillCommand()                 
        if( jobId instanceof Collection ) {             
            result.add(jobId.join(','))
    //        log.trace "Kill command: ${result}"         
        }                                               
        else {                                          
            result.add(jobId.toString())                
        }                                               
        return result                                   
    }                                                   

    @Override
    protected List<String> queueStatusCommand(Object queue) {
        return TupleHelper.listOf('hq','--output-mode=json','job', 'list', '--all')
    }

    @Override
    protected Map<String, QueueStatus> parseQueueStatus(String text) {
        final jobs = parser.parseText(text) as List<Map>
        final result = new HashMap()
        for( Map entry : jobs ) {
            final id = entry.id as String
            final status = parseJobStatus(entry)
            result.put(id, status)
        }
        return result
    }

    protected QueueStatus parseJobStatus( Map entry ) {
        Map stats = entry.task_stats as Map
        if( !stats )
            return QueueStatus.UNKNOWN
        if( stats.canceled )
            return QueueStatus.DONE
        if( stats.failed )
            return QueueStatus.ERROR
        if( stats.finished )
            return QueueStatus.DONE
        if( stats.running )
            return QueueStatus.RUNNING
        if( stats.waiting )
            return QueueStatus.HOLD
        else
            return QueueStatus.UNKNOWN
    }
}
