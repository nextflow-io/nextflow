/*
 * Copyright 2022, Google Inc.
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

package nextflow.cloud.google.batch

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.cloud.google.batch.client.BatchClient
import nextflow.cloud.google.batch.client.BatchConfig
import nextflow.cloud.google.batch.model.AllocationPolicy
import nextflow.cloud.google.batch.model.BatchJob
import nextflow.cloud.google.batch.model.ComputeResource
import nextflow.cloud.google.batch.model.NetworkInterface
import nextflow.cloud.google.batch.model.NetworkPolicy
import nextflow.cloud.google.batch.model.ProvisioningModel
import nextflow.cloud.google.batch.model.TaskContainer
import nextflow.cloud.google.batch.model.TaskGroup
import nextflow.cloud.google.batch.model.TaskRunnable
import nextflow.cloud.google.batch.model.TaskSpec
import nextflow.processor.TaskBean
import nextflow.processor.TaskHandler
import nextflow.processor.TaskRun
import nextflow.processor.TaskStatus
/**
 * Implements a task handler for Google Batch executor
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class GoogleBatchTaskHandler extends TaskHandler {

    private GoogleBatchExecutor executor

    private TaskBean taskBean

    private Path exitFile

    private Path outputFile

    private Path errorFile

    private BatchClient client

    private String jobId

    private String jobState

    private long timestamp

    private GoogleBatchScriptLauncher launcher

    GoogleBatchTaskHandler(TaskRun task, GoogleBatchExecutor executor) {
        super(task)
        this.client = executor.getClient()
        this.jobId = "nf-${task.hashLog.replace('/','')}-${System.currentTimeMillis()}"
        this.executor = executor
        this.taskBean = task.toTaskBean()
        this.launcher = new GoogleBatchScriptLauncher(taskBean, executor.remoteBinDir)
        // those files are access via NF runtime, keep based on CloudStoragePath
        this.outputFile = task.workDir.resolve(TaskRun.CMD_OUTFILE)
        this.errorFile = task.workDir.resolve(TaskRun.CMD_ERRFILE)
        this.exitFile = task.workDir.resolve(TaskRun.CMD_EXIT)
    }


    /*
     * Only for testing -- do not use
     */
    protected GoogleBatchTaskHandler() {}

    @Override
    void submit() {
        /*
         * create the task runner script
         */
        launcher.build()

        /*
         * create submit request
         */
        final req = newSubmitRequest(task)
        log.debug "[GOOGLE BATCH] new job request > $req"
        final resp = client.submitJob(jobId, req)
        this.status = TaskStatus.SUBMITTED
        log.debug "[GOOGLE BATCH] submitted > job=$jobId; work-dir=${task.getWorkDirStr()}; resp=$resp"
    }

    protected BatchJob newSubmitRequest(TaskRun task) {
        final result = new BatchJob()

        final spec = new TaskSpec()
        final res = new ComputeResource()
        // CPUs requirement
        res.cpuMilli = task.config.getCpus() * 1000
        // memory requirement
        if( task.config.getMemory() )
            res.memoryMib = task.config.getMemory().getMega().toInteger()
        // timeout requirement
        if( task.config.getTime() )
            spec.withMaxRunDuration(task.config.getTime() )

        // task spec
        final cmd = "trap \"{ cp ${TaskRun.CMD_LOG} ${launcher.workDirMount}/${TaskRun.CMD_LOG}; }\" ERR; /bin/bash ${launcher.workDirMount}/${TaskRun.CMD_RUN} 2>&1 | tee ${TaskRun.CMD_LOG}"
        final container = new TaskContainer()
                .withImageUri(task.container)
                .withCommands(['/bin/bash','-o','pipefail','-c', cmd.toString()])
                // note: container must mount the base work directory
                .withVolumes(launcher.getContainerMounts())
                .withOptions(task.config.getContainerOptions())

        spec.addRunnable(new TaskRunnable(container: container))
                .withVolumes(launcher.getTaskVolumes())
                .withComputeResources(res)

        // allocation policy
        final allocPolicy = new AllocationPolicy()
        final networkPolicy= networkPolicy(executor.config)
        if( networkPolicy )
            allocPolicy.withNetworkPolicy(networkPolicy)
        if( task.config.getMachineType() )
            allocPolicy.withMachineTypes( task.config.getMachineType() )
        if( executor.config.preemptible )
            allocPolicy.withProvisioningModel(ProvisioningModel.PREEMPTIBLE)
        if( executor.config.spot )
            allocPolicy.withProvisioningModel(ProvisioningModel.SPOT)
        
        // create the task group
        result.addTaskGroup( new TaskGroup().withTaskSpec(spec).withAllocationPolicy(allocPolicy) )
        return result
    }

    protected NetworkPolicy networkPolicy(BatchConfig config) {
        NetworkInterface result=null
        // set the network
        if( config.network ) {
            if(!result) result = new NetworkInterface()
            result.withNetwork(config.network)
        }
        // set the subnetwork
        if( config.subnetwork ) {
            if(!result) result = new NetworkInterface()
            result.withSubnetwork(config.subnetwork)
        }
        // set private address
        if( config.usePrivateAddress ) {
            if(!result) result = new NetworkInterface()
            result.withNoExternalIpAddress(true)
        }

        return result
                ? new NetworkPolicy().withNetworkInterfaces(result)
                : null
    }
    /**
     * @return Retrieve the submitted pod state
     */
    protected String getJobState() {
        final now = System.currentTimeMillis()
        final delta =  now - timestamp;
        if( !jobState || delta >= 1_000) {
            def newState = client.getJobState(jobId)
            if( newState ) {
                log.trace "[GOOGLE BATCH] Get Batch job=$jobId state=$newState"
                jobState = newState
                timestamp = now
            }
        }
        return jobState
    }

    private List<String> RUNNING_AND_TERMINATED = ['RUNNING', 'SUCCEEDED', 'FAILED', 'DELETION_IN_PROGRESS']

    private List<String> TERMINATED = ['SUCCEEDED', 'FAILED', 'DELETION_IN_PROGRESS']


    @Override
    boolean checkIfRunning() {
        if(isSubmitted()) {
            // include `terminated` state to allow the handler status to progress
            if (getJobState() in RUNNING_AND_TERMINATED) {
                status = TaskStatus.RUNNING
                return true
            }
        }
        return false
    }

    @Override
    boolean checkIfCompleted() {
        if( getJobState() in TERMINATED ) {
            // finalize the task
            task.exitStatus = readExitFile()
            task.stdout = outputFile
            task.stderr = errorFile
            status = TaskStatus.COMPLETED
            return true
        }

        return false
    }

    @PackageScope Integer readExitFile() {
        try {
            exitFile.text as Integer
        }
        catch (Exception e) {
            log.debug "[GOOGLE BATCH] Cannot read exitstatus for task: `$task.name` | ${e.message}"
            null
        }
    }

    @Override
    void kill() {
        if( isSubmitted() ) {
            log.trace "[GOOGLE BATCH] deleting pod name=$jobId"
            client.deleteJobs(jobId)
        }
        else {
            log.debug "[GOOGLE BATCH] Oops.. invalid delete action"
        }
    }

}
