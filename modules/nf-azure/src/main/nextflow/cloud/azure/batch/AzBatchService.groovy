/*
 * Copyright 2020, Microsoft Corp
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

package nextflow.cloud.azure.batch

import com.microsoft.azure.batch.BatchClient
import com.microsoft.azure.batch.auth.BatchSharedKeyCredentials
import com.microsoft.azure.batch.protocol.models.CloudTask
import com.microsoft.azure.batch.protocol.models.ContainerConfiguration
import com.microsoft.azure.batch.protocol.models.ImageInformation
import com.microsoft.azure.batch.protocol.models.OutputFile
import com.microsoft.azure.batch.protocol.models.PoolInformation
import com.microsoft.azure.batch.protocol.models.ResourceFile
import com.microsoft.azure.batch.protocol.models.TaskAddParameter
import com.microsoft.azure.batch.protocol.models.TaskContainerSettings
import com.microsoft.azure.batch.protocol.models.VirtualMachineConfiguration
import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.util.logging.Slf4j
import nextflow.cloud.azure.config.AzConfig
import nextflow.cloud.azure.config.AzPoolOpts
import nextflow.cloud.azure.file.AzStorageContainerParser
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
/**
 * Implements Azure Batch operations for Nextflow executor
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class AzBatchService {

    static final private List<String> allPoolIds = new ArrayList<>(50)

    AzPoolOpts poolOpts

    AzConfig config

    Map<TaskProcessor,String> allJobIds = new HashMap<>(50)

    AzBatchService(AzBatchExecutor executor) {
        assert executor
        this.config = executor.config
        this.poolOpts = config.batch().pool()
    }

    @Memoized
    protected BatchClient getClient() {
        createBatchClient()
    }

    protected BatchClient createBatchClient() {
        log.debug "[AZURE BATCH] Executor options=${config.batch()}"
        // Create batch client
        if( !config.batch().endpoint )
            throw new IllegalArgumentException("Missing Azure Batch endpoint -- Specify it in the nextflow.config file using the setting 'azure.batch.endpoint'")
        if( !config.batch().accountName )
            throw new IllegalArgumentException("Missing Azure Batch account name -- Specify it in the nextflow.config file using the setting 'azure.batch.accountName'")
        if( !config.batch().accountKey )
            throw new IllegalArgumentException("Missing Azure Batch account key -- Specify it in the nextflow.config file using the setting 'azure.batch.accountKet'")

        final cred = new BatchSharedKeyCredentials(config.batch().endpoint, config.batch().accountName, config.batch().accountKey)
        return BatchClient.open(cred)
    }


    AzTaskKey submitTask(TaskRun task) {
        final poolId = getOrCreatePool()
        final jobId = getOrCreateJob(poolId, task)
        runTask(jobId, task)
    }

    CloudTask getTask(AzTaskKey key) {
        return client.taskOperations().getTask(key.jobId, key.taskId)
    }

    void terminate(AzTaskKey key) {
        client.taskOperations().terminateTask(key.jobId, key.taskId)
    }

    synchronized String getOrCreateJob(String poolId, TaskRun task) {
        if( allJobIds.containsKey(task.processor)) {
            return allJobIds[task]
        }
        // create a batch job
        final jobId = makeJobId(task)
        final poolInfo = new PoolInformation()
        poolInfo.withPoolId(poolId)
        client.jobOperations().createJob(jobId, poolInfo)
        // add to the map
        allJobIds[task.processor] = jobId
        return jobId
    }

    String makeJobId(TaskRun task) {
        def name = task
                .processor
                .name.trim().replaceAll(/[^a-zA-Z0-9-_]+/,'_')
        return "nf-job-$name"
    }

    AzTaskKey runTask(String jobId, TaskRun task) {
        final container = task.config.container as String
        if( !container )
            throw new IllegalArgumentException("Missing container image for process: $task.name")

        final taskId = "nf-${task.hash.toString()}"
        final containerOpts = new TaskContainerSettings().withImageName(container)

        TaskAddParameter taskToAdd = new TaskAddParameter()
                .withId(taskId)
                .withContainerSettings(containerOpts)
                .withCommandLine("bash ${TaskRun.CMD_RUN}")
                .withResourceFiles(resourceFiles(task))
                .withOutputFiles(outputFiles(task))
        client.taskOperations().createTask(jobId, taskToAdd)
        return new AzTaskKey(jobId, taskId)
    }

    protected List<ResourceFile> resourceFiles(TaskRun task) {
        final uri = task.workDir.toUriString()
        final bucket = AzStorageContainerParser.parse(uri)
        if( !bucket )
            throw new IllegalArgumentException("Invalid Azure data container URI: $uri")

        final result = new ArrayList(10)

        // add the bash launcher + script
        final launcher = new ResourceFile()
                .withAutoStorageContainerName(bucket.container)
                .withBlobPrefix("$bucket.key/${TaskRun.CMD_RUN}")
                .withFilePath("$bucket.path/${TaskRun.CMD_RUN}")
                .withFileMode('554')
        result.add(launcher)

        final script = new ResourceFile()
                .withAutoStorageContainerName(bucket.container)
                .withBlobPrefix("$bucket.key/${TaskRun.CMD_SCRIPT}")
                .withFilePath("$bucket.path/${TaskRun.CMD_SCRIPT}")
                .withFileMode('554')
        result.add(script)

        return result
    }

    protected List<OutputFile> outputFiles(TaskRun task) {
        return []
    }

    protected ImageInformation getImage() {
        List<ImageInformation> images = client.accountOperations().listSupportedImages()

        for (ImageInformation it : images) {
            if( it.osType() != poolOpts.type )
                continue
            if( it.verificationType() != poolOpts.verification )
                continue
            if( !it.imageReference().publisher().equalsIgnoreCase(poolOpts.publisher) )
                continue
            if( it.imageReference().offer().equalsIgnoreCase(poolOpts.offer) )
                return it
        }

        throw new IllegalStateException("Cannot find a matching VM image with publister=$poolOpts.publisher; offer=$poolOpts.offer; OS type=$poolOpts.type; verification type=$poolOpts.verification")
    }

    synchronized String getOrCreatePool() {
        // define a stable pool name
        final String poolId = "nf-pool-$poolOpts.vmType"

        // check existence and create if needed
        if( !client.poolOperations().existsPool(poolId) ) {

            final image = getImage()

            /**
             * A container configuration must be provided for a task to run in a specific container.
             * Such container can be pre-fetched on VM creation or when running the task
             *
             * https://github.com/MicrosoftDocs/azure-docs/blob/master/articles/batch/batch-docker-container-workloads.md#:~:text=Run%20container%20applications%20on%20Azure,compatible%20containers%20on%20the%20nodes.
             */
            final containerConfig = new ContainerConfiguration();

            final vmConfig = new VirtualMachineConfiguration()
                    .withNodeAgentSKUId(image.nodeAgentSKUId())
                    .withImageReference(image.imageReference())
                    .withContainerConfiguration(containerConfig)

            client.poolOperations().createPool(poolId, poolOpts.vmType, vmConfig, poolOpts.vmCount)
            // add to the list of pool ids
            allPoolIds.add(poolId)
        }

        return poolId
    }



}
