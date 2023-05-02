/*
 * Copyright 2013-2023, Seqera Labs
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

package nextflow.cloud.aws.batch

import static nextflow.fusion.FusionConfig.FUSION_PATH
import static nextflow.fusion.FusionHelper.*

import com.amazonaws.services.batch.AWSBatch
import com.amazonaws.services.batch.model.ArrayProperties
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.executor.TaskArraySubmitter
import nextflow.processor.TaskHandler
import nextflow.processor.TaskRun
import nextflow.processor.TaskStatus
/**
 * Submit tasks as an array job for a grid executor.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class AwsBatchTaskArraySubmitter extends TaskArraySubmitter implements SubmitJobAware {

    private AwsBatchExecutor executor

    private AWSBatch client

    AwsBatchTaskArraySubmitter(List<TaskHandler> array, AwsBatchExecutor executor) {
        super(array)
        this.executor = executor
        this.client = executor.client

        if( array.size() > 10_000 )
            throw new IllegalArgumentException("Array jobs on AWS Batch may not have size greater than 10,000")
    }

    @Override
    TaskRun getTask() { array.first().getTask() }

    @Override
    AWSBatch getClient() { executor.client }

    @Override
    AwsOptions getAwsOptions() { executor.getAwsOptions() }

    @Override
    List<String> getSubmitCommand() {
        // create wrapper script
        final arrayIndexName = executor.getArrayIndexName()
        final workDirs = fusionEnabled()
            ? array.collect { h -> toContainerMount(h.task.workDir).toString() }
            : array.collect { h -> h.task.workDir.toUriString() }

        def cmd = ((AwsBatchTaskHandler)array.first())
            .getSubmitCommand()
            .last()
            .replaceAll(workDirs.first(), '\\$task_dir')

        if( fusionEnabled() )
            cmd = "bash ${cmd}"

        cmd = """
            declare -a array=( ${workDirs.join(' ')} )
            task_dir=\${array[\$${arrayIndexName}]}
            ${cmd}
            """.stripIndent().trim()

        // create command line
        final cli = ['bash', '-o', 'pipefail', '-c', cmd.toString()]
        if( fusionEnabled() )
            cli.add(0, FUSION_PATH)

        return cli
    }

    @Override
    void submit() {
        // -- create the submit request
        final request = newSubmitRequest(getTask())
            .withArrayProperties(new ArrayProperties().withSize(array.size()))

        // -- submit the array job
        final response = submitJobRequest(bypassProxy(client), request)
        final jobId = response.jobId

        // -- set the job id, queue, and status of each task
        array.eachWithIndex { handler, i ->
            ((AwsBatchTaskHandler)handler).jobId = executor.getArrayTaskId(jobId, i)
            ((AwsBatchTaskHandler)handler).queueName = request.getJobQueue()
            handler.status = TaskStatus.SUBMITTED
        }

        log.debug "[AWS BATCH] submitted array job > jobId: ${jobId}"
    }

}
