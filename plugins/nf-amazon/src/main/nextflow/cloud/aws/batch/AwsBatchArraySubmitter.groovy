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
import com.amazonaws.services.batch.model.AWSBatchException
import com.amazonaws.services.batch.model.ArrayProperties
import com.amazonaws.services.batch.model.SubmitJobRequest
import com.amazonaws.services.batch.model.SubmitJobResult
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.exception.ProcessSubmitException
import nextflow.executor.ArrayTaskSubmitter
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
class AwsBatchArraySubmitter extends ArrayTaskSubmitter implements SubmitJobAware {

    private AwsBatchExecutor executor

    private AWSBatch client

    AwsBatchArraySubmitter(List<TaskHandler> array, AwsBatchExecutor executor) {
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
        return fusionEnabled()
                ? fusionSubmitCli()
                : classicSubmitCli()
    }

    @Override
    List<String> fusionSubmitCli() {
        final workDirs = array
            .collect { handler -> toContainerMount(handler.task.workDir) }
            .join(' ')

        final cmd = """
            declare -a array=( ${workDirs} )
            bash \${array[\$AWS_BATCH_JOB_ARRAY_INDEX]}/${TaskRun.CMD_RUN}
            """.stripIndent().trim()

        return List.of(FUSION_PATH, 'bash', cmd.toString())
    }

    protected List<String> classicSubmitCli() {
        final workDirs = array
            .collect { handler -> handler.task.workDir.toUriString() }
            .join(' ')

        final opts = getAwsOptions()
        final cli = opts.getAwsCli()
        final debug = opts.debug ? ' --debug' : ''
        final sse = opts.storageEncryption ? " --sse ${opts.storageEncryption}" : ''
        final kms = opts.storageKmsKeyId ? " --sse-kms-key-id ${opts.storageKmsKeyId}" : ''
        final aws = "${cli} s3 cp --only-show-errors${debug}${sse}${kms}"

        final cmd = """
            declare -a array=( ${workDirs} )
            task_dir=\${array[\$AWS_BATCH_JOB_ARRAY_INDEX]}
            trap "{ ret=\$?; ${aws} ${TaskRun.CMD_LOG} \$task_dir/${TaskRun.CMD_LOG}||true; exit \$ret; }" EXIT
            ${aws} \$task_dir/${TaskRun.CMD_RUN} - | bash 2>&1 | tee ${TaskRun.CMD_LOG}
            """.stripIndent().trim()

        return List.of('bash', '-o', 'pipefail', '-c', cmd.toString())
    }

    @Override
    void submit() {
        // -- create the submit request
        final request = newSubmitRequest(getTask())
            .withArrayProperties(new ArrayProperties().withSize(array.size()))

        // -- submit the array job
        final response = submitJobRequest0(bypassProxy(client), request)
        final jobId = response.jobId

        // -- set the job id, queue, and status of each task
        array.eachWithIndex { handler, i ->
            ((AwsBatchTaskHandler)handler).jobId = getArrayTaskId(jobId, i)
            ((AwsBatchTaskHandler)handler).queueName = request.getJobQueue()
            handler.status = TaskStatus.SUBMITTED
        }

        log.debug "[AWS BATCH] submitted array job > jobId: ${jobId}"
    }

    static private SubmitJobResult submitJobRequest0(AWSBatch client, SubmitJobRequest request) {
        try {
            return client.submitJob(request)
        }
        catch( AWSBatchException e ) {
            if( e.statusCode >= 500 )
                // raise a process exception so that nextflow can try to recover it
                throw new ProcessSubmitException("Failed to submit job: ${request.jobName} - Reason: ${e.errorCode}", e)
            else
                // status code < 500 are not expected to be recoverable, just throw it again
                throw e
        }
    }
 
    String getArrayTaskId(String jobId, int index) {
        "${jobId}:${index}"
    }

}
