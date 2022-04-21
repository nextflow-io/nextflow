/*
 * Copyright (c) 2020-2021. Seqera Labs, S.L.
 *
 * All Rights reserved
 *
 */

package nextflow.cloud.google.batch

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.cloud.google.batch.client.BatchClient
import nextflow.cloud.google.batch.model.BatchJob
import nextflow.cloud.google.batch.model.ComputeResource
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
        this.taskBean = new TaskBean(task)
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
        // create the task group
        result.addTaskGroup( new TaskGroup().withTaskSpec(spec) )
        return result
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
//            savePodLogOnError(task)
//            deletePodIfSuccessful(task)
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
            log.trace "[K8s] deleting pod name=$jobId"
            client.deleteJobs(jobId)
        }
        else {
            log.debug "[K8s] Oops.. invalid delete action"
        }
    }

}
