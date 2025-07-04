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
import io.seqera.client.SeqeraClient
import io.seqera.sched.api.v1a1.CreateJobRequest
import io.seqera.sched.api.v1a1.GetJobLogsResponse
import io.seqera.sched.api.v1a1.JobState
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

    private SeqeraClient client

    private SeqeraExecutor executor

    private Path exitFile

    private Path outputFile

    private Path errorFile

    private String jobId

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
        final req = new CreateJobRequest()
            .withCommand(fusionSubmitCli())
            .withImage(task.getContainer())
            .withClusterId(executor.getClusterId())
            .withEnvironment(fusionLauncher().fusionEnv())
            .withPlatform(task.getContainerPlatform())
            .withCpus(task.config.getCpus())
            .withMemory(task.config.getMemory().toBytes())
        log.debug "[SEQERA] Submitting job request=${req}"
        final resp = client.createJob(req)
        this.jobId = resp.jobId
        this.status = TaskStatus.SUBMITTED
    }

    protected JobState jobState() {
        return client
            .describeJob(jobId)
            .jobsState
    }

    @Override
    boolean checkIfRunning() {
        if (isSubmitted()) {
            final state = jobState()
            log.debug "[SEQERA] checkIfRunning job=${jobId}; state=${state}"
            if (state.isRunningOrTerminated()) {
                status = TaskStatus.RUNNING
                return true
            }
        }
        return false
    }

    @Override
    boolean checkIfCompleted() {
        final state = jobState()
        log.debug "[SEQERA] checkIfCompleted state=${state}"
        if (state.isTerminated()) {
            log.debug "[SEQERA] Process `${task.lazyName()}` - terminated job=$jobId; state=$state"
            // finalize the task
            task.exitStatus = readExitFile()
            if (state.isFailed()) {
                final logs = getJobLogs(jobId)
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

    protected GetJobLogsResponse getJobLogs(String jobId) {
        return client.getJobLogs(jobId)
    }

    @Override
    protected void killTask() {
        log.debug "[SEQERA] Kill job=${jobId}"
        client.cancelJob(jobId)
    }

    @PackageScope
    Integer readExitFile() {
        try {
            final result = exitFile.text as Integer
            log.trace "[SEQERA] Read exit file for job $jobId; exit=${result}"
            return result
        }
        catch (Exception e) {
            log.debug "[SEQERA] Cannot read exit status for task: `${task.lazyName()}` - ${e.message}"
            // return MAX_VALUE to signal it was unable to retrieve the exit code
            return Integer.MAX_VALUE
        }
    }
}
