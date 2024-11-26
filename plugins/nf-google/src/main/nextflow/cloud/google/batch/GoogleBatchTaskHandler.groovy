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
import java.util.regex.Pattern

import com.google.cloud.batch.v1.AllocationPolicy
import com.google.cloud.batch.v1.ComputeResource
import com.google.cloud.batch.v1.Environment
import com.google.cloud.batch.v1.Job
import com.google.cloud.batch.v1.LifecyclePolicy
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
import nextflow.exception.ProcessException
import nextflow.exception.ProcessUnrecoverableException
import nextflow.executor.BashWrapperBuilder
import nextflow.executor.res.DiskResource
import nextflow.fusion.FusionAwareTask
import nextflow.fusion.FusionScriptLauncher
import nextflow.processor.TaskArrayRun
import nextflow.processor.TaskConfig
import nextflow.processor.TaskHandler
import nextflow.processor.TaskProcessor
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

    private static Pattern EXIT_CODE_REGEX = ~/exit code 500(\d\d)/

    private GoogleBatchExecutor executor

    private Path exitFile

    private Path outputFile

    private Path errorFile

    private BatchClient client

    private BashWrapperBuilder launcher

    /**
     * Job Id assigned by Nextflow
     */
    private String jobId

    /**
     * Task id assigned by Google Batch service
     */
    private String taskId

    /**
     * Job unique id assigned by Google Batch service
     */
    private String uid

    /**
     * Task state assigned by Google Batch service
     */
    private String taskState

    private volatile CloudMachineInfo machineInfo

    private volatile long timestamp

    GoogleBatchTaskHandler(TaskRun task, GoogleBatchExecutor executor) {
        super(task)
        this.client = executor.getClient()
        this.executor = executor
        this.jobId = customJobName(task) ?: "nf-${task.hashLog.replace('/','')}-${System.currentTimeMillis()}"
        // those files are access via NF runtime, keep based on CloudStoragePath
        this.outputFile = task.workDir.resolve(TaskRun.CMD_OUTFILE)
        this.errorFile = task.workDir.resolve(TaskRun.CMD_ERRFILE)
        this.exitFile = task.workDir.resolve(TaskRun.CMD_EXIT)
    }

    /**
     * Resolve the `jobName` property defined in the nextflow config file
     *
     * @param task The underlying task to be executed
     * @return The custom job name for the specified task or {@code null} if the `jobName` attribute has been specified
     */
    protected String customJobName(TaskRun task) {
        try {
            final custom = (Closure)executor.session?.getExecConfigProp(executor.name, 'jobName', null)
            if( !custom )
                return null

            final ctx = [ (TaskProcessor.TASK_CONTEXT_PROPERTY_NAME): task.config ]
            return custom.cloneWith(ctx).call()?.toString()
        }
        catch( Exception e ) {
            log.debug "Unable to resolve job custom name", e
            return null
        }
    }

    protected BashWrapperBuilder createTaskWrapper() {
        if( fusionEnabled() ) {
            return fusionLauncher()
        }
        else {
            final taskBean = task.toTaskBean()
            return new GoogleBatchScriptLauncher(taskBean, executor.remoteBinDir)
                .withConfig(executor.config)
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
    void prepareLauncher() {
        launcher = createTaskWrapper()
        launcher.build()
    }

    @Override
    void submit() {
        /*
         * create submit request
         */
        final req = newSubmitRequest(task, spec0(launcher))
        log.trace "[GOOGLE BATCH] new job request > $req"
        final resp = client.submitJob(jobId, req)
        final uid = resp.getUid()
        updateStatus(jobId, '0', uid)
        log.debug "[GOOGLE BATCH] Process `${task.lazyName()}` submitted > job=$jobId; uid=$uid; work-dir=${task.getWorkDirStr()}"
    }

    protected void updateStatus(String jobId, String taskId, String uid) {
        if( task instanceof TaskArrayRun ) {
            // update status for children
            for( int i=0; i<task.children.size(); i++ ) {
                final handler = task.children[i] as GoogleBatchTaskHandler
                final arrayTaskId = executor.getArrayTaskId(jobId, i)
                handler.updateStatus(jobId, arrayTaskId, uid)
            }
        }
        else {
            this.jobId = jobId
            this.taskId = taskId
            this.uid = uid
            this.status = TaskStatus.SUBMITTED
        }
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
        if( accel && accel.type.toLowerCase().startsWith('nvidia-') ) {
            container
                .addVolumes('/var/lib/nvidia/lib64:/usr/local/nvidia/lib64')
                .addVolumes('/var/lib/nvidia/bin:/usr/local/nvidia/bin')
        }

        def containerOptions = task.config.getContainerOptions() ?: ''
        // accelerator requires privileged option
        // https://cloud.google.com/batch/docs/create-run-job#create-job-gpu
        if( task.config.getAccelerator() || fusionEnabled() ) {
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

        // retry on spot reclaim
        if( executor.config.maxSpotAttempts ) {
            // Note: Google Batch uses the special exit status 50001 to signal
            // the execution was terminated due a spot reclaim. When this happens
            // The policy re-execute the jobs automatically up to `maxSpotAttempts` times
            taskSpec
                .setMaxRetryCount( executor.config.maxSpotAttempts )
                .addLifecyclePolicies(
                    LifecyclePolicy.newBuilder()
                        .setActionCondition(
                            LifecyclePolicy.ActionCondition.newBuilder()
                                .addExitCodes(50001)
                        )
                        .setAction(LifecyclePolicy.Action.RETRY_TASK)
                )
        }

        // instance policy
        // allocation policy
        final allocationPolicy = AllocationPolicy.newBuilder()
        final instancePolicyOrTemplate = AllocationPolicy.InstancePolicyOrTemplate.newBuilder()

        if( executor.config.getAllowedLocations() )
            allocationPolicy.setLocation(
                AllocationPolicy.LocationPolicy.newBuilder()
                    .addAllAllowedLocations( executor.config.getAllowedLocations() )
            )

        if( executor.config.serviceAccountEmail )
            allocationPolicy.setServiceAccount(
                ServiceAccount.newBuilder()
                    .setEmail( executor.config.serviceAccountEmail )
            )

        allocationPolicy.putAllLabels( task.config.getResourceLabels() )

        // use instance template if specified
        if( task.config.getMachineType()?.startsWith('template://') ) {
            if( task.config.getAccelerator() )
                log.warn1 'Process directive `accelerator` ignored because an instance template was specified'

            if( task.config.getDisk() )
                log.warn1 'Process directive `disk` ignored because an instance template was specified'

            if( executor.config.getBootDiskImage() )
                log.warn1 'Config option `google.batch.bootDiskImage` ignored because an instance template was specified'

            if( executor.config.cpuPlatform )
                log.warn1 'Config option `google.batch.cpuPlatform` ignored because an instance template was specified'

            if( executor.config.preemptible )
                log.warn1 'Config option `google.batch.premptible` ignored because an instance template was specified'

            if( executor.config.spot )
                log.warn1 'Config option `google.batch.spot` ignored because an instance template was specified'

            instancePolicyOrTemplate
                .setInstallGpuDrivers( executor.config.getInstallGpuDrivers() )
                .setInstanceTemplate( task.config.getMachineType().minus('template://') )
        }

        // otherwise create instance policy
        else {
            final instancePolicy = AllocationPolicy.InstancePolicy.newBuilder()

            if( task.config.getAccelerator() ) {
                final accelerator = AllocationPolicy.Accelerator.newBuilder()
                    .setCount( task.config.getAccelerator().getRequest() )

                if( task.config.getAccelerator().getType() )
                    accelerator.setType( task.config.getAccelerator().getType() )

                instancePolicy.addAccelerators(accelerator)
                instancePolicyOrTemplate.setInstallGpuDrivers(true)
            }

            if( executor.config.getBootDiskImage() )
                instancePolicy.setBootDisk( AllocationPolicy.Disk.newBuilder().setImage( executor.config.getBootDiskImage() ) )

            if( fusionEnabled() && !disk ) {
                disk = new DiskResource(request: '375 GB', type: 'local-ssd')
                log.debug "[GOOGLE BATCH] Process `${task.lazyName()}` - adding local volume as fusion scratch: $disk"
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

            if( executor.config.cpuPlatform )
                instancePolicy.setMinCpuPlatform( executor.config.cpuPlatform )

            if( executor.config.preemptible )
                instancePolicy.setProvisioningModel( AllocationPolicy.ProvisioningModel.PREEMPTIBLE )

            if( executor.config.spot )
                instancePolicy.setProvisioningModel( AllocationPolicy.ProvisioningModel.SPOT )

            instancePolicyOrTemplate.setPolicy( instancePolicy )
        }

        allocationPolicy.addInstances(instancePolicyOrTemplate)

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

        // task group
        final taskGroup = TaskGroup.newBuilder()
            .setTaskSpec(taskSpec)

        if( task instanceof TaskArrayRun ) {
            final arraySize = task.getArraySize()
            taskGroup.setTaskCount(arraySize)
        }

        // create the job
        return Job.newBuilder()
            .addTaskGroups(taskGroup)
            .setAllocationPolicy(allocationPolicy)
            .setLogsPolicy(
                LogsPolicy.newBuilder()
                    .setDestination(LogsPolicy.Destination.CLOUD_LOGGING)
            )
            .putAllLabels(task.config.getResourceLabels())
            .build()
    }

    /**
     * @return Retrieve the submitted task state
     */
    protected String getTaskState() {
        final tasks = client.listTasks(jobId)
        if( !tasks.iterator().hasNext() )
            return 'PENDING'

        final now = System.currentTimeMillis()
        final delta =  now - timestamp;
        if( !taskState || delta >= 1_000) {
            final status = client.getTaskStatus(jobId, taskId)
            final newState = status?.state as String
            if( newState ) {
                log.trace "[GOOGLE BATCH] Get job=$jobId task=$taskId state=$newState"
                taskState = newState
                timestamp = now
            }
            if( newState == 'PENDING' ) {
                final eventsCount = status.getStatusEventsCount()
                final lastEvent = eventsCount > 0 ? status.getStatusEvents(eventsCount - 1) : null
                if( lastEvent?.getDescription()?.contains('CODE_GCE_QUOTA_EXCEEDED') )
                    log.warn1 "Batch job cannot be run: ${lastEvent.getDescription()}"
            }
        }
        return taskState
    }

    static private final List<String> RUNNING_OR_COMPLETED = ['RUNNING', 'SUCCEEDED', 'FAILED']

    static private final List<String> COMPLETED = ['SUCCEEDED', 'FAILED']

    @Override
    boolean checkIfRunning() {
        if(isSubmitted()) {
            // include `terminated` state to allow the handler status to progress
            if( getTaskState() in RUNNING_OR_COMPLETED ) {
                status = TaskStatus.RUNNING
                return true
            }
        }
        return false
    }

    @Override
    boolean checkIfCompleted() {
        final state = getTaskState()
        if( state in COMPLETED ) {
            log.debug "[GOOGLE BATCH] Process `${task.lazyName()}` - terminated job=$jobId; task=$taskId; state=$state"
            // finalize the task
            task.exitStatus = readExitFile()
            if( state == 'FAILED' ) {
                if( task.exitStatus == Integer.MAX_VALUE )
                    task.error = getJobError()
                task.stdout = executor.logging.stdout(uid, taskId) ?: outputFile
                task.stderr = executor.logging.stderr(uid, taskId) ?: errorFile
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

    protected Throwable getJobError() {
        try {
            final status = client.getTaskStatus(jobId, taskId)
            final eventsCount = status.getStatusEventsCount()
            final lastEvent = eventsCount > 0 ? status.getStatusEvents(eventsCount - 1) : null
            log.debug "[GOOGLE BATCH] Process `${task.lazyName()}` - last event: ${lastEvent}; exit code: ${lastEvent?.taskExecution?.exitCode}"

            final error = lastEvent?.description
            if( error && EXIT_CODE_REGEX.matcher(error).find() ) {
                return new ProcessException(error)
            }
        }
        catch (Throwable t) {
            log.debug "[GOOGLE BATCH] Unable to fetch task `${task.lazyName()}` exit code - cause: ${t.message}"
        }

        return null
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
        if( isActive() ) {
            log.trace "[GOOGLE BATCH] Process `${task.lazyName()}` - deleting job name=$jobId"
            if( executor.shouldDeleteJob(jobId) )
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
        if( jobId && uid )
            result.put('native_id', "$jobId/$taskId/$uid")
        result.machineInfo = getMachineInfo()
        return result
    }

    protected GoogleBatchMachineTypeSelector.MachineType bestMachineType0(int cpus, int memory, String location, boolean spot, boolean localSSD, List<String> families) {
        return GoogleBatchMachineTypeSelector.INSTANCE.bestMachineType(cpus, memory, location, spot, localSSD, families)
    }

    protected GoogleBatchMachineTypeSelector.MachineType findBestMachineType(TaskConfig config, boolean localSSD) {
        final location = client.location
        final cpus = config.getCpus()
        final memory = config.getMemory() ? config.getMemory().toMega().toInteger() : 1024
        final spot = executor.config.spot ?: executor.config.preemptible
        final machineType = config.getMachineType()
        final families = machineType ? machineType.tokenize(',') : List.<String>of()
        final priceModel = spot ? PriceModel.spot : PriceModel.standard

        try {
            if( executor.isCloudinfoEnabled() ) {
                return bestMachineType0(cpus, memory, location, spot, localSSD, families)
            }
        }
        catch (Exception e) {
            log.debug "[GOOGLE BATCH] Cannot select machine type using Seqera Cloudinfo for task: `${task.lazyName()}` - ${e.message}"
        }

        // Check if a specific machine type was provided by the user
        if( machineType && !machineType.contains(',') && !machineType.contains('*') )
            return new GoogleBatchMachineTypeSelector.MachineType(
                type: machineType,
                location: location,
                priceModel: priceModel
            )

        // Fallback to Google Batch automatically deduce from requested resources
        return null
    }

}
