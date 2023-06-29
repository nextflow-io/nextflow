/*
 * Copyright 2023, Seqera Labs
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

import com.google.cloud.batch.v1.AllocationPolicy
import com.google.cloud.batch.v1.ComputeResource
import com.google.cloud.batch.v1.Environment
import com.google.cloud.batch.v1.Job
import com.google.cloud.batch.v1.LogsPolicy
import com.google.cloud.batch.v1.Runnable
import com.google.cloud.batch.v1.ServiceAccount
import com.google.cloud.batch.v1.TaskGroup
import com.google.cloud.batch.v1.TaskSpec
import com.google.cloud.batch.v1.Volume
import com.google.protobuf.Duration
import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.cloud.google.batch.client.BatchClient
import nextflow.cloud.types.CloudMachineInfo
import nextflow.cloud.types.PriceModel
import nextflow.exception.ProcessUnrecoverableException
import nextflow.executor.BashWrapperBuilder
import nextflow.executor.res.DiskResource
import nextflow.fusion.FusionAwareTask
import nextflow.fusion.FusionScriptLauncher
import nextflow.processor.TaskConfig
import nextflow.processor.TaskHandler
import nextflow.processor.TaskRun
import nextflow.processor.TaskStatus
import nextflow.trace.TraceRecord
/**
 * Implements a task handler for Google Batch executor
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class GoogleBatchTaskHandler extends TaskHandler implements FusionAwareTask {

    private GoogleBatchExecutor executor

    private Path exitFile

    private Path outputFile

    private Path errorFile

    private BatchClient client

    /**
     * Job Id assigned by Nextflow
     */
    private String jobId

    /**
     * Job unique id assigned by Google Batch service
     */
    private String uid

    /**
     * Job state assigned by Google Batch service
     */
    private String jobState

    private volatile CloudMachineInfo machineInfo

    private volatile long timestamp

    GoogleBatchTaskHandler(TaskRun task, GoogleBatchExecutor executor) {
        super(task)
        this.client = executor.getClient()
        this.jobId = "nf-${task.hashLog.replace('/','')}-${System.currentTimeMillis()}"
        this.executor = executor
        // those files are access via NF runtime, keep based on CloudStoragePath
        this.outputFile = task.workDir.resolve(TaskRun.CMD_OUTFILE)
        this.errorFile = task.workDir.resolve(TaskRun.CMD_ERRFILE)
        this.exitFile = task.workDir.resolve(TaskRun.CMD_EXIT)
    }

    protected BashWrapperBuilder createTaskWrapper() {
        if( fusionEnabled() ) {
            return fusionLauncher()
        }
        else {
            final taskBean = task.toTaskBean()
            return new GoogleBatchScriptLauncher(taskBean, executor.remoteBinDir)
        }
    }

    /*
     * Only for testing -- do not use
     */
    protected GoogleBatchTaskHandler() {}

    protected GoogleBatchLauncherSpec spec0(BashWrapperBuilder launcher) {
        if( launcher instanceof GoogleBatchLauncherSpec )
            return launcher
        if( launcher instanceof FusionScriptLauncher )
            return new GoogleBatchFusionAdapter(this, launcher)
        throw new IllegalArgumentException("Unexpected Google Batch launcher type: ${launcher?.getClass()?.getName()}")
    }

    @Override
    void submit() {
        /*
         * create the task runner script
         */
        final launcher = createTaskWrapper()
        launcher.build()

        /*
         * create submit request
         */
        final req = newSubmitRequest(task, spec0(launcher))
        log.trace "[GOOGLE BATCH] new job request > $req"
        final resp = client.submitJob(jobId, req)
        this.uid = resp.getUid()
        this.status = TaskStatus.SUBMITTED
        log.debug "[GOOGLE BATCH] Process `${task.lazyName()}` submitted > job=$jobId; uid=$uid; work-dir=${task.getWorkDirStr()}"
    }

    protected Job newSubmitRequest(TaskRun task, GoogleBatchLauncherSpec launcher) {
        // resource requirements
        final taskSpec = TaskSpec.newBuilder()
        final computeResource = ComputeResource.newBuilder()

        computeResource.setCpuMilli( task.config.getCpus() * 1000 )

        if( task.config.getMemory() )
            computeResource.setMemoryMib( task.config.getMemory().getMega() )

        if( task.config.getTime() )
            taskSpec.setMaxRunDuration(
                Duration.newBuilder()
                    .setSeconds( task.config.getTime().toSeconds() )
            )

        def disk = task.config.getDiskResource()
        // apply disk directive to boot disk if type is not specified
        if( disk && !disk.type )
            computeResource.setBootDiskMib( disk.request.getMega() )
        // otherwise use config setting
        else if( executor.config.bootDiskSize )
            computeResource.setBootDiskMib( executor.config.bootDiskSize.getMega() )

        // container
        if( !task.container )
            throw new ProcessUnrecoverableException("Process `${task.lazyName()}` failed because the container image was not specified")

        final cmd = launcher.launchCommand()
        final container = Runnable.Container.newBuilder()
            .setImageUri( task.container )
            .addAllCommands( cmd )
            .addAllVolumes( launcher.getContainerMounts() )

        final accel = task.config.getAccelerator()
        // add nvidia specific driver paths
        // see https://cloud.google.com/batch/docs/create-run-job#create-job-gpu
        if(  accel && accel.type.toLowerCase().startsWith('nvidia-') ) {
            container
                .addVolumes('/var/lib/nvidia/lib64:/usr/local/nvidia/lib64')
                .addVolumes('/var/lib/nvidia/bin:/usr/local/nvidia/bin')
        }

        def containerOptions= task.config.getContainerOptions() ?: ''
        // accelerator requires privileged option
        // https://cloud.google.com/batch/docs/create-run-job#create-job-gpu
        if( task.config.getAccelerator() || fusionEnabled()) {
            if( containerOptions ) containerOptions += ' '
            containerOptions += '--privileged'
        }

        if( containerOptions )
            container.setOptions( containerOptions )

        // task spec
        final env = Environment
                .newBuilder()
                .putAllVariables( launcher.getEnvironment() )
                .build()

        taskSpec
            .setComputeResource(computeResource)
            .addRunnables(
                Runnable.newBuilder()
                    .setContainer(container)
                    .setEnvironment(env)
            )
            .addAllVolumes( launcher.getVolumes() )

        // instance policy
        final allocationPolicy = AllocationPolicy.newBuilder()
        final instancePolicyOrTemplate = AllocationPolicy.InstancePolicyOrTemplate.newBuilder()
        final instancePolicy = AllocationPolicy.InstancePolicy.newBuilder()

        if( executor.config.getAllowedLocations() )
            allocationPolicy.setLocation(
                AllocationPolicy.LocationPolicy.newBuilder()
                    .addAllAllowedLocations( executor.config.getAllowedLocations() )
            )

        if( task.config.getAccelerator() ) {
            final accelerator = AllocationPolicy.Accelerator.newBuilder()
                .setCount( task.config.getAccelerator().getRequest() )

            if( task.config.getAccelerator().getType() )
                accelerator.setType( task.config.getAccelerator().getType() )

            instancePolicy.addAccelerators(accelerator)
            instancePolicyOrTemplate.setInstallGpuDrivers(true)
        }

        if( fusionEnabled() && !disk ) {
            disk = new DiskResource(request: '375 GB', type: 'local-ssd')
            log.debug "[GOOGLE BATCH] Process `${task.lazyName()}` - adding local volume as fusion scratch: $disk"
        }

        if( executor.config.cpuPlatform ) {
            instancePolicy.setMinCpuPlatform( executor.config.cpuPlatform )
        }

        final machineType = findBestMachineType(task.config, disk?.type == 'local-ssd')
        if( machineType ) {
            instancePolicy.setMachineType(machineType.type)
            machineInfo = new CloudMachineInfo(
                    type: machineType.type,
                    zone: machineType.location,
                    priceModel: machineType.priceModel
            )
        }

        // When using local SSD not all the disk sizes are valid and depends on the machine type
        if( disk?.type == 'local-ssd' && machineType ) {
            final validSize = GoogleBatchMachineTypeSelector.INSTANCE.findValidLocalSSDSize(disk.request, machineType)
            if( validSize != disk.request ) {
                disk = new DiskResource(request: validSize, type: 'local-ssd')
                log.debug "[GOOGLE BATCH] Process `${task.lazyName()}` - adjusting local disk size to: $validSize"
            }
        }

        // use disk directive for an attached disk if type is specified
        if( disk?.type ) {
            instancePolicy.addDisks(
                AllocationPolicy.AttachedDisk.newBuilder()
                    .setNewDisk(
                        AllocationPolicy.Disk.newBuilder()
                            .setType(disk.type)
                            .setSizeGb(disk.request.toGiga())
                    )
                    .setDeviceName('scratch')
            )

            taskSpec.addVolumes(
                Volume.newBuilder()
                    .setDeviceName('scratch')
                    .setMountPath('/tmp')
            )
        }

        if( executor.config.serviceAccountEmail )
            allocationPolicy.setServiceAccount(
                ServiceAccount.newBuilder()
                    .setEmail(executor.config.serviceAccountEmail)
            )

        if( executor.config.preemptible )
            instancePolicy.setProvisioningModel( AllocationPolicy.ProvisioningModel.PREEMPTIBLE )

        if( executor.config.spot )
            instancePolicy.setProvisioningModel( AllocationPolicy.ProvisioningModel.SPOT )

        allocationPolicy.addInstances(
            instancePolicyOrTemplate
                .setPolicy(instancePolicy)
        )

        // network policy
        final networkInterface = AllocationPolicy.NetworkInterface.newBuilder()
        def hasNetworkPolicy = false

        if( executor.config.network ) {
            hasNetworkPolicy = true
            networkInterface.setNetwork( executor.config.network )
        }
        if( executor.config.subnetwork ) {
            hasNetworkPolicy = true
            networkInterface.setSubnetwork( executor.config.subnetwork )
        }
        if( executor.config.usePrivateAddress ) {
            hasNetworkPolicy = true
            networkInterface.setNoExternalIpAddress( true )
        }

        if( hasNetworkPolicy )
            allocationPolicy.setNetwork(
                AllocationPolicy.NetworkPolicy.newBuilder()
                    .addNetworkInterfaces(networkInterface)
            )

        allocationPolicy.putAllLabels(task.config.getResourceLabels())

        // create the job
        return Job.newBuilder()
            .addTaskGroups(
                TaskGroup.newBuilder()
                    .setTaskSpec(taskSpec)
            )
            .setAllocationPolicy(allocationPolicy)
            .setLogsPolicy(
                LogsPolicy.newBuilder()
                    .setDestination(LogsPolicy.Destination.CLOUD_LOGGING)
            )
            .build()
    }

    /**
     * @return Retrieve the submitted job state
     */
    protected String getJobState() {
        final now = System.currentTimeMillis()
        final delta =  now - timestamp;
        if( !jobState || delta >= 1_000) {
            final status = client.getJobStatus(jobId)
            final newState = status?.state as String
            if( newState ) {
                log.trace "[GOOGLE BATCH] Get job=$jobId state=$newState"
                jobState = newState
                timestamp = now
            }
            if( newState == 'SCHEDULED' ) {
                final eventsCount = status.getStatusEventsCount()
                final lastEvent = eventsCount > 0 ? status.getStatusEvents(eventsCount - 1) : null
                if( lastEvent?.getDescription()?.contains('CODE_GCE_QUOTA_EXCEEDED') )
                    log.warn1 "Batch job cannot be run: ${lastEvent.getDescription()}"
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
        final state = getJobState()
        if( state in TERMINATED ) {
            log.debug "[GOOGLE BATCH] Process `${task.lazyName()}` - terminated job=$jobId; state=$state"
            // finalize the task
            task.exitStatus = readExitFile()
            if( state == 'FAILED' ) {
                task.stdout = executor.logging.stdout(uid) ?: outputFile
                task.stderr = executor.logging.stderr(uid) ?: errorFile
            }
            else {
                task.stdout = outputFile
                task.stderr = errorFile
            }
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
            log.debug "[GOOGLE BATCH] Cannot read exit status for task: `${task.lazyName()}` - ${e.message}"
            // return MAX_VALUE to signal it was unable to retrieve the exit code
            return Integer.MAX_VALUE
        }
    }

    @Override
    void kill() {
        if( isSubmitted() ) {
            log.trace "[GOOGLE BATCH] Process `${task.lazyName()}` - deleting job name=$jobId"
            client.deleteJob(jobId)
        }
        else {
            log.debug "[GOOGLE BATCH] Process `${task.lazyName()}` - invalid delete action"
        }
    }

    protected CloudMachineInfo getMachineInfo() {
        return machineInfo
    }

    @Override
    TraceRecord getTraceRecord() {
        def result = super.getTraceRecord()
        if( jobId && uid ) {
            result.put('native_id', "$jobId/$uid")
        }
        result.machineInfo = getMachineInfo()
        return result
    }

    protected GoogleBatchMachineTypeSelector.MachineType findBestMachineType(TaskConfig config, boolean localSSD) {
        final location = client.location
        final cpus = config.getCpus()
        final memory = config.getMemory() ? config.getMemory().toMega().toInteger() : 1024
        final spot = executor.config.spot ?: executor.config.preemptible
        final families = config.getMachineType() ? config.getMachineType().tokenize(',') : []
        final priceModel = spot ? PriceModel.spot : PriceModel.standard

        try {
            return GoogleBatchMachineTypeSelector.INSTANCE.bestMachineType(cpus, memory, location, spot, localSSD, families)
        }
        catch (Exception e) {
            log.debug "[GOOGLE BATCH] Cannot select machine type using cloud info for task: `${task.lazyName()}` - ${e.message}"

            // Check if a specific machine type was provided by the user
            if( config.getMachineType() && !config.getMachineType().contains(',') && !config.getMachineType().contains('*') )
                return new GoogleBatchMachineTypeSelector.MachineType(
                        type: config.getMachineType(),
                        location: location,
                        priceModel: priceModel
                )

            // Fallback to Google Batch automatically deduce from requested resources
            return null
        }

    }

}
