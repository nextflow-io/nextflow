/*
 * Copyright 2013-2025, Seqera Labs
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

package io.seqera.executor

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import io.seqera.sched.api.schema.v1a1.AcceleratorType
import io.seqera.sched.api.schema.v1a1.GetTaskLogsResponse
import io.seqera.sched.api.schema.v1a1.NextflowTask
import io.seqera.sched.api.schema.v1a1.ResourceLimit
import io.seqera.sched.api.schema.v1a1.ResourceRequirement
import io.seqera.sched.api.schema.v1a1.Task
import io.seqera.sched.api.schema.v1a1.TaskState as SchedTaskState
import io.seqera.sched.api.schema.v1a1.TaskStatus as SchedTaskStatus
import io.seqera.sched.client.SchedClient
import io.seqera.util.MapperUtil
import nextflow.cloud.types.CloudMachineInfo
import nextflow.exception.ProcessException
import nextflow.exception.ProcessUnrecoverableException
import nextflow.util.Duration
import nextflow.util.MemoryUnit
import nextflow.fusion.FusionAwareTask
import nextflow.processor.TaskHandler
import nextflow.processor.TaskRun
import nextflow.processor.TaskStatus
import nextflow.trace.TraceRecord
/**
 * Task handler for the Seqera scheduler executor.
 *
 * <p>Manages the lifecycle of a single task submitted to the Seqera scheduler,
 * including submission via batch submitter, status polling, completion handling,
 * and trace record enrichment with machine info and spot interruption metadata.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class SeqeraTaskHandler extends TaskHandler implements FusionAwareTask {

    private SchedClient client

    private SeqeraExecutor executor

    private Path exitFile

    private Path outputFile

    private Path errorFile

    private volatile String taskId

    /**
     * Cached task state from last describeTask call, used for trace record metadata
     */
    private volatile SchedTaskState cachedTaskState

    /**
     * Cached machine info extracted from task attempts
     */
    private volatile CloudMachineInfo machineInfo

    SeqeraTaskHandler(TaskRun task, SeqeraExecutor executor) {
        super(task)
        this.client = executor.getClient()
        this.executor = executor
        // those files are access via NF runtime, keep based on CloudStoragePath
        this.outputFile = task.workDir.resolve(TaskRun.CMD_OUTFILE)
        this.errorFile = task.workDir.resolve(TaskRun.CMD_ERRFILE)
        this.exitFile = task.workDir.resolve(TaskRun.CMD_EXIT)
    }

    @Override
    void prepareLauncher() {
        assert fusionEnabled()
        final launcher = fusionLauncher()
        launcher.build()
    }

    @Override
    void submit() {
        executor.ensureRunCreated()
        int cpuShares = (task.config.getCpus() ?: 1) * 1024
        int memoryMiB = task.config.getMemory() ? (int) (task.config.getMemory().toBytes() / (1024 * 1024)) : 1024
        final resourceReq = new ResourceRequirement()
            .cpuShares(cpuShares)
            .memoryMiB(memoryMiB)
        // add accelerator settings if defined
        final accelerator = task.config.getAccelerator()
        if( accelerator ) {
            // number of accelerators requested, fallback to limit if request is not specified
            resourceReq.acceleratorCount(accelerator.request ?: accelerator.limit)
            // accelerator type is GPU by default (most common in scientific computing)
            resourceReq.acceleratorType(AcceleratorType.GPU)
            // specific accelerator model name e.g. "nvidia-tesla-v100", "nvidia-a10g"
            if( accelerator.type )
                resourceReq.acceleratorName(accelerator.type)
        }
        // build machine requirement merging config settings with task arch, disk, and snapshot settings
        final machineReq = MapperUtil.toMachineRequirement(
            executor.getSeqeraConfig().machineRequirement,
            task.getContainerPlatform(),
            task.config.getDisk(),
            fusionConfig().snapshotsEnabled()
        )
        // build resource limit from process resourceLimits directive (upper bound for OOM retry scaling)
        final resourceLim = toResourceLimit()
        // validate container - Seqera executor requires all processes to specify a container image
        final container = task.getContainer()
        if( !container )
            throw new ProcessUnrecoverableException("Process `${task.lazyName()}` failed because the container image was not specified -- the Seqera executor requires all processes define a container image")
        // build the scheduler task with all required attributes
        final schedTask = new Task()
            .name(task.lazyName())       // process name for identification
            .image(container)             // container image to run
            .command(fusionSubmitCli())   // fusion-based command launcher
            .environment(fusionLauncher().fusionEnv())  // fusion environment variables
            .resourceRequirement(resourceReq)  // cpu, memory, accelerators
            .resourceLimit(resourceLim)         // resource upper bounds for OOM retry
            .machineRequirement(machineReq)    // machine type and disk requirements
            .nextflow(new NextflowTask()
                .taskId(task.id?.intValue())
                .hash(task.hash?.toString())
                .workDir(task.getWorkDirStr()))
        log.debug "[SEQERA] Enqueueing task for batch submission: ${schedTask}"
        // Enqueue for batch submission - status will be set by setBatchTaskId callback
        executor.getBatchSubmitter().submit(this, schedTask)
    }

    /**
     * Called by batch submitter after successful batch submission
     */
    void setBatchTaskId(String taskId) {
        this.taskId = taskId
        this.status = TaskStatus.SUBMITTED
        log.debug "[SEQERA] Process `${task.lazyName()}` submitted > taskId=$taskId; work-dir=${task.getWorkDirStr()}"
    }

    /**
     * Called by batch submitter when batch submission fails
     */
    void onBatchSubmitFailure(Exception cause) {
        log.debug "[SEQERA] Batch submission failed for task ${task.lazyName()}: ${cause.message}"
        task.error = cause
        this.status = TaskStatus.COMPLETED
    }

    /**
     * Build a {@link ResourceLimit} from the process {@code resourceLimits} directive.
     * Returns {@code null} if no resource limits are defined.
     */
    protected ResourceLimit toResourceLimit() {
        final memoryLimit = task.config.getResourceLimit('memory') as MemoryUnit
        final cpusLimit = task.config.getResourceLimit('cpus') as Integer
        if( !memoryLimit && !cpusLimit )
            return null
        final result = new ResourceLimit()
        if( memoryLimit )
            result.memoryMiB((int)(memoryLimit.toBytes() / (1024 * 1024)))
        if( cpusLimit )
            result.cpuShares(cpusLimit * 1024)
        return result
    }

    protected SchedTaskStatus schedTaskStatus() {
        cachedTaskState = client.describeTask(taskId).getTaskState()
        return cachedTaskState.getStatus()
    }

    @Override
    boolean checkIfRunning() {
        if (isSubmitted()) {
            final schedStatus = schedTaskStatus()
            log.debug "[SEQERA] checkIfRunning taskId=${taskId}; status=${schedStatus}"
            if (isRunningOrTerminated(schedStatus)) {
                status = TaskStatus.RUNNING
                return true
            }
        }
        return false
    }

    @Override
    boolean checkIfCompleted() {
        // Handle batch submission failure - task error was set but never reached RUNNING state
        if (task.error && isCompleted()) {
            return true
        }
        if (!isRunning())
            return false
        final schedStatus = schedTaskStatus()
        log.debug "[SEQERA] checkIfCompleted status=${schedStatus}"
        if (isTerminated(schedStatus)) {
            log.debug "[SEQERA] Process `${task.lazyName()}` - terminated taskId=$taskId; status=$schedStatus"
            // finalize the task
            task.exitStatus = readExitFile()
            if (isFailed(schedStatus)) {
                // When no exit code available, get the error message from task state
                if (task.exitStatus == Integer.MAX_VALUE) {
                    final errorMessage = cachedTaskState?.getErrorMessage() ?: "Task failed for unknown reason"
                    task.error = new ProcessException(errorMessage)
                }
                final logs = getTaskLogs(taskId)
                task.stdout = logs?.stdout ?: outputFile
                task.stderr = logs?.stderr ?: errorFile
            } else {
                task.stdout = outputFile
                task.stderr = errorFile
            }
            status = TaskStatus.COMPLETED
            return true
        }

        return false
    }

    protected boolean isRunningOrTerminated(SchedTaskStatus status) {
        return status == SchedTaskStatus.RUNNING || isTerminated(status)
    }

    protected boolean isTerminated(SchedTaskStatus status) {
        return status in [SchedTaskStatus.SUCCEEDED, SchedTaskStatus.FAILED, SchedTaskStatus.CANCELLED]
    }

    protected boolean isFailed(SchedTaskStatus status) {
        return status == SchedTaskStatus.FAILED
    }

    protected GetTaskLogsResponse getTaskLogs(String taskId) {
        return client.getTaskLogs(taskId)
    }

    @Override
    protected void killTask() {
        if( !taskId ) {
            log.trace "[SEQERA] Skip cancel - taskId not yet assigned"
            return
        }
        log.debug "[SEQERA] Cancel taskId=${taskId}"
        try {
            client.cancelTask(taskId)
        }
        catch (Throwable t) {
            log.warn "[SEQERA] Failed to cancel task ${taskId}", t
        }
    }

    @PackageScope
    Integer readExitFile() {
        try {
            final result = exitFile.text as Integer
            log.trace "[SEQERA] Read exit file for taskId $taskId; exit=${result}"
            return result
        }
        catch (Exception e) {
            log.debug "[SEQERA] Cannot read exit status for task: `${task.lazyName()}` - ${e.message}"
            // return MAX_VALUE to signal it was unable to retrieve the exit code
            return Integer.MAX_VALUE
        }
    }

    /**
     * Get machine info for the task execution from the last task attempt.
     * The machine info is cached after first retrieval.
     *
     * @return CloudMachineInfo containing instance type, zone, and price model, or null if not available
     */
    protected CloudMachineInfo getMachineInfo() {
        if (machineInfo)
            return machineInfo
        if (!cachedTaskState)
            return null

        try {
            final attempts = cachedTaskState.getAttempts()
            if (!attempts || attempts.isEmpty())
                return null

            final lastAttempt = attempts.get(attempts.size() - 1)
            final lastInfo = lastAttempt.getMachineInfo()
            if (!lastInfo)
                return null

            // Convert Sched API MachineInfo to Nextflow CloudMachineInfo
            machineInfo = new CloudMachineInfo(
                type: lastInfo.getType(),
                zone: lastInfo.getZone(),
                priceModel: MapperUtil.toPriceModel(lastInfo.getPriceModel())
            )
            log.trace "[SEQERA] taskId=$taskId => machineInfo=$machineInfo"
            return machineInfo
        }
        catch (Exception e) {
            log.debug "[SEQERA] Unable to get machine info for taskId=$taskId - ${e.message}"
            return null
        }
    }

    /**
     * Get the number of spot interruptions for this task.
     * This is calculated server-side from task attempts with spot-related stop reasons.
     *
     * @return the count of spot interruptions, or null if not completed or not available
     */
    protected Integer getNumSpotInterruptions() {
        if (!taskId || !isCompleted())
            return null
        if (!cachedTaskState)
            return null
        return cachedTaskState.getNumSpotInterruptions()
    }

    /**
     * Get the native backend ID for this task (ECS task ARN or Docker container ID).
     *
     * @return the native ID from the last task attempt, or null if not available
     */
    protected String getNativeId() {
        return cachedTaskState?.getId()
    }

    protected Long getGrantedTime() {
        final time = cachedTaskState?.getResourceRequirement()?.getTime()
        return time != null ? Duration.of(time).toMillis() : task.config.getTime()?.toMillis()
    }

    /**
     * Get the trace record for this task, including machine info and spot interruptions metadata.
     *
     * @return the trace record with additional metadata fields
     */
    @Override
    TraceRecord getTraceRecord() {
        final result = super.getTraceRecord()
        result.put('native_id', getNativeId())
        result.machineInfo = getMachineInfo()
        result.numSpotInterruptions = getNumSpotInterruptions()
        // Override executor name to include cloud backend for cost tracking
        result.executorName = "${SeqeraExecutor.SEQERA}/aws"
        return result
    }
}
