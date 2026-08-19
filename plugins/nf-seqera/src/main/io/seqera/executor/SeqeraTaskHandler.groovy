/*
 * Copyright 2013-2026, Seqera Labs
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

package io.seqera.executor

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import io.seqera.config.MachineRequirementOpts
import io.seqera.executor.Labels
import io.seqera.sched.api.schema.v1a1.AcceleratorType
import io.seqera.sched.api.schema.v1a1.GetTaskLogsResponse
import io.seqera.sched.api.schema.v1a1.NextflowTask
import io.seqera.sched.api.schema.v1a1.PredictionModel
import io.seqera.sched.api.schema.v1a1.ResourceLimit
import io.seqera.sched.api.schema.v1a1.ResourceRequirement
import io.seqera.sched.api.schema.v1a1.Task
import io.seqera.sched.api.schema.v1a1.TaskState as SchedTaskState
import io.seqera.sched.api.schema.v1a1.TaskStatus as SchedTaskStatus
import io.seqera.sched.client.SchedClient
import io.seqera.util.HintHelper
import io.seqera.util.SchemaMapperUtil
import nextflow.cloud.types.CloudMachineInfo
import nextflow.exception.ProcessException
import nextflow.exception.ProcessUnrecoverableException
import nextflow.util.Duration
import nextflow.util.MemoryUnit
import nextflow.fusion.FusionAwareTask
import nextflow.fusion.FusionConfig
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
        // cpus needs no unspecified handling: TaskConfig.getCpus() already guarantees at least 1
        // and throws on a negative, so there is never an absent or non-positive value to forward.
        // memory has no such guarantee — see memoryMiB().
        final resourceReq = new ResourceRequirement()
            .cpuShares(task.config.getCpus() * 1024)
            .memoryMiB(memoryMiB())
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
        // overlay any seqera/machineRequirement.* hints on top of config-scope values (hints win)
        final baseMachineOpts = HintHelper.overlayHints(
            executor.getSeqeraConfig().machineRequirement,
            task.config.getHints()
        )
        final machineReq = SchemaMapperUtil.toMachineRequirement(
            baseMachineOpts,
            task.getContainerPlatform(),
            task.config.getDisk(),
            fusionConfig().snapshotsEnabled(),
            maxSpotAttempts(baseMachineOpts)
        )
        // Resolve the optional per-task prediction model override from the seqera/predictionModel
        // hint. An explicit hint always wins: the automatic check below is a safety net, and the
        // user asking for a specific model on a process is a deliberate opt-out of it.
        // When neither applies the value is left null, so the task inherits the run-level model
        final predictionModelHint = HintHelper.resolvePredictionModel(task.config.getHints())
        final predictionModel = predictionModelHint
            ? PredictionModel.fromValue(predictionModelHint)
            : (shouldDisablePrediction() ? PredictionModel.NONE : null)
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
            .environment(getTaskEnvironment())  // fusion + user-configured environment variables
            .resourceRequirement(resourceReq)  // cpu, memory, accelerators
            .resourceLimit(resourceLim)         // resource upper bounds for OOM retry
            .machineRequirement(machineReq)    // machine type and disk requirements
            .predictionModel(predictionModel)  // optional per-task prediction model override
            .nextflow(new NextflowTask()
                .taskId(task.id?.intValue())
                .hash(task.hash?.toString())
                .workDir(task.getWorkDirStr()))
        // attach per-task resource labels delta (over run-level baseline)
        final taskLabels = Labels.toStringMap(task.config.getResourceLabels())
        final delta = Labels.delta(taskLabels, executor.runResourceLabels)
        if( delta )
            schedTask.labels(delta)
        // attach pipeline secret references (never values) — resolved at the compute edge
        final secretRefs = buildSecretRefs()
        if( secretRefs )
            schedTask.secrets(secretRefs)
        log.debug "[SEQERA] Enqueueing task for batch submission: ${schedTask}"
        // Enqueue for batch submission - status will be set by setBatchTaskId callback
        executor.getBatchSubmitter().submit(this, schedTask)
    }

    /**
     * Determine whether the resource prediction model must be disabled for this task.
     *
     * The prediction model can allocate less memory than the task requested. However the task script
     * is rendered *before* the task is scheduled, therefore a script referencing {@code task.memory}
     * carries the memory that was requested, not the one that has been allocated e.g. a JVM
     * {@code -Xmx} setting exceeding the container memory and failing with an out-of-memory error.
     *
     * @return {@code true} when the task should be submitted with prediction model {@code none}
     */
    protected boolean shouldDisablePrediction() {
        // Nothing to disable unless the run enables a prediction model. Checking this first
        // also keeps the warning quiet for the runs where the resources are never adjusted,
        // which would otherwise report a problem that cannot happen
        final runModel = executor.getSeqeraConfig()?.predictionModel
        if( !runModel || runModel == PredictionModel.NONE.getValue() )
            return false
        // Note only `memory` is checked. The scheduler can adjust the cpus as well, but a
        // stale `task.cpus` costs an over-subscribed thread pool whereas a stale
        // `task.memory` fails the task, and `task.cpus` is referenced by nearly every
        // process -- disabling the prediction for those would defeat the feature
        if( !task.isDirectiveReferenced('memory') )
            return false
        // warn once per process rather than once per task: the reference belongs to the
        // process definition, so every one of its tasks would otherwise report it
        log.warn1("Process `${task.processor?.name ?: task.lazyName()}` depends on the `task.memory` value -- resource prediction has been disabled for this process to prevent an under-allocation of the requested memory", firstOnly: true)
        return true
    }

    protected int maxSpotAttempts(MachineRequirementOpts opts) {
        final result = opts?.maxSpotAttempts
        if( result != null && result < 0 )
            throw new IllegalArgumentException("Invalid maxSpotAttempts value: ${result} -- the value must be zero or a positive number")
        if( result )
            return result
        // when fusion snapshot is enabled max attempt should be > 0
        // to enable to allow snapshot retry the job execution in a new compute instance
        return fusionConfig().snapshotsEnabled() ? FusionConfig.DEFAULT_SNAPSHOT_MAX_SPOT_ATTEMPTS : 0
    }

    /**
     * Build the map of container environment variable name to the pipeline secret store
     * reference for each {@code secret} process directive.
     *
     * Platform stores the secret value in the cloud secret store under
     * {@code tower-<workflowId>/<name>}; the executor sends only this reference and the
     * scheduler backend resolves the value at the compute edge. The secret value never
     * passes through Nextflow or the scheduler API.
     *
     * Returns {@code null} when there are no secret directives or no Platform workflow id
     * (bare-scheduler usage has no store, so a missing reference is a no-op).
     */
    protected Map<String, String> buildSecretRefs() {
        final names = task.config.getSecret()
        final workflowId = executor.getWorkflowId()
        if( !names || !workflowId )
            return null
        final result = new LinkedHashMap<String, String>()
        for( String name : names )
            result.put(name, "tower-${workflowId}/${name}".toString())
        return result
    }

    /**
     * Build the task environment by merging user-configured environment variables
     * with Fusion environment variables. Fusion variables take precedence.
     */
    protected Map<String, String> getTaskEnvironment() {
        final configEnv = executor.getSeqeraConfig()?.taskEnvironment
        final fusionEnv = fusionLauncher().fusionEnv()
        if( !configEnv )
            return fusionEnv
        final result = new LinkedHashMap<String, String>(configEnv)
        result.putAll(fusionEnv)
        return result
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
     * The memory to request, in MiB, or {@code null} when the process declares none.
     * <p>
     * The {@code memory} directive is optional, and leaving it out is a legitimate way of saying
     * "I have no opinion". Substituting a figure here fabricated an opinion the user never
     * expressed: the scheduler received a concrete request indistinguishable from a declared one,
     * so it could neither apply its own default nor report that nothing had been asked for, and a
     * whole run of unspecified tasks looked deliberately sized. Omitting the field instead lets the
     * scheduler resolve and surface its default (seqeralabs/sched#1086).
     * <p>
     * A zero-valued directive is omitted for the same reason, and warned about once per process,
     * since unlike an absent directive it is a config accident rather than a choice. Exactly zero
     * is the only invalid figure reachable here: {@code MemoryUnit(long)} asserts a non-negative
     * value, and {@code TaskConfig.getMemory0} maps a falsy directive to null via
     * {@code MemoryUnit.asBoolean}, so what survives is the string form — {@code memory '0 GB'} —
     * whose string is truthy.
     * <p>
     * The conversion below can still yield zero for a <em>positive</em> sub-MiB directive: the
     * size-derived idiom nf-core uses for index builds, {@code memory { 6.B * fasta.size() }} in
     * {@code nf-core/bwa/index}, falls under 1 MiB for any reference below ~170 KB — which is what
     * a test or CI profile supplies, and is why {@code BWAMEM1_INDEX} submits a zero today. That
     * zero is forwarded deliberately rather than caught here. The scheduler already resolves any
     * non-positive request to its default and reports having done so; restating that rule
     * client-side would split one normalisation across two codebases that then have to be kept in
     * agreement, which is the failure this whole change is undoing.
     *
     * @return the requested memory in MiB, or null to leave the axis unspecified
     */
    protected Integer memoryMiB() {
        final memory = task.config.getMemory()
        if( memory == null )
            return null
        if( memory.toBytes() <= 0 ) {
            log.warn1("Process `${task.processor.name}` declares a zero `memory` directive -- ignoring it; the scheduler will apply its own default", firstOnly: true)
            return null
        }
        return (int) (memory.toBytes() / (1024 * 1024))
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
            // prefer the exit code reported by the scheduler API; fall back to the `.exitcode`
            // file only when the API does not report one. On error (e.g. OOM, spot reclaim,
            // timeout) the container may terminate before the wrapper's on_exit trap can write
            // the file, so the scheduler exit code is the more reliable source — consistent with
            // the K8s, AWS Batch and Azure Batch executors.
            final apiExitCode = cachedTaskState?.getExitCode()
            task.exitStatus = apiExitCode != null ? apiExitCode : readExitFile()
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
                priceModel: SchemaMapperUtil.toPriceModel(lastInfo.getPriceModel())
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
     * Get the log stream identifier for this task.
     *
     * @return the log stream ID, or null if not available
     */
    protected String getLogStreamId() {
        return cachedTaskState?.getLogStreamId()
    }

    /**
     * Get the native backend ID for this task (ECS task ARN or Docker container ID).
     *
     * @return the native ID from the last task attempt, or null if not available
     */
    protected String getNativeId() {
        return cachedTaskState?.getId()
    }

    /**
     * Get the allocated resources for this task from the last task attempt.
     * Falls back to the resource requirement from the task state if no attempts exist.
     *
     * @return a map of allocated resource fields, or null if not available
     */
    protected Map<String,Object> getResourceAllocation() {
        if (!cachedTaskState)
            return null

        def resources = null
        final attempts = cachedTaskState.getAttempts()
        if (attempts && !attempts.isEmpty()) {
            resources = attempts.get(attempts.size() - 1).getResources()
        }
        if (!resources) {
            resources = cachedTaskState.getResourceAllocation()
        }
        if (!resources)
            return null

        final result = new LinkedHashMap<String,Object>()
        if (resources.getCpuShares() != null)
            result.put('cpuShares', resources.getCpuShares())
        if (resources.getMemoryMiB() != null)
            result.put('memoryMiB', resources.getMemoryMiB())
        if (resources.getAcceleratorCount() != null)
            result.put('acceleratorCount', resources.getAcceleratorCount())
        if (resources.getAcceleratorType() != null)
            result.put('acceleratorType', resources.getAcceleratorType().toString())
        if (resources.getAcceleratorName() != null)
            result.put('acceleratorName', resources.getAcceleratorName())
        if (resources.getTime() != null)
            result.put('time', resources.getTime())
        return result.isEmpty() ? null : result
    }

    protected Long getGrantedTime() {
        final String time = cachedTaskState?.getResourceAllocation()?.getTime()
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
        result.logStreamId = getLogStreamId()
        result.resourceAllocation = getResourceAllocation()
        // Override executor name to include cloud backend for cost tracking
        result.executorName = "${SeqeraExecutor.SEQERA}/aws"
        return result
    }
}
