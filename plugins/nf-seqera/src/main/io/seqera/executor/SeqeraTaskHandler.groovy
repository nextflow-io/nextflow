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

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import io.seqera.sched.api.schema.v1a1.GetTaskLogsResponse
import io.seqera.sched.api.schema.v1a1.Task
import io.seqera.sched.api.schema.v1a1.TaskStatus as SchedTaskStatus
import io.seqera.sched.client.SchedClient
import nextflow.fusion.FusionAwareTask
import nextflow.processor.TaskHandler
import nextflow.processor.TaskRun
import nextflow.processor.TaskStatus

import java.nio.file.Path

/**
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

    private String taskId

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
        int cpuRequest = task.config.getCpus() ?: 1
        long memoryRequest = task.config.getMemory() ? task.config.getMemory().toBytes() : 1024 * 1024 * 1024
        final schedTask = new Task()
            .sessionId(executor.getSessionId())
            .image(task.getContainer())
            .command(fusionSubmitCli())
            .environment(fusionLauncher().fusionEnv())
            .arch(task.getContainerPlatform())
            .cpus(cpuRequest)
            .memory(memoryRequest)
        log.debug "[SEQERA] Enqueueing task for batch submission: ${schedTask}"
        // Enqueue for batch submission - status will be set by setBatchTaskId callback
        executor.getBatchSubmitter().enqueue(this, schedTask)
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

    protected SchedTaskStatus schedTaskStatus() {
        return client
            .describeTask(taskId)
            .getTaskState()
            .getStatus()
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
        final schedStatus = schedTaskStatus()
        log.debug "[SEQERA] checkIfCompleted status=${schedStatus}"
        if (isTerminated(schedStatus)) {
            log.debug "[SEQERA] Process `${task.lazyName()}` - terminated taskId=$taskId; status=$schedStatus"
            // finalize the task
            task.exitStatus = readExitFile()
            if (isFailed(schedStatus)) {
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
        log.debug "[SEQERA] Cancel taskId=${taskId}"
        client.cancelTask(taskId)
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
}
