/*
 * Copyright 2019, Genome Research Limited
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

package nextflow.wr.executor

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import java.nio.file.Path

import nextflow.wr.client.WrRestApi

import nextflow.processor.TaskRun
import nextflow.processor.TaskHandler
import nextflow.processor.TaskStatus
import nextflow.processor.BatchContext
import nextflow.processor.BatchHandler

/**
 * Handles a job execution using wr as a backend, without needing files on disk
 * for state.
 *
 * @author Sendu Bala <sb10@sanger.ac.uk>
 * Based on TesTaskHandler by Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class WrTaskHandler extends TaskHandler implements BatchHandler<String,Map> {

    public static final List<String> COMPLETE_STATUSES = ['complete', 'buried', 'deleted']
    public static final List<String> STARTED_STATUSES = ['delayed', 'reserved', 'running'] + COMPLETE_STATUSES

    final WrExecutor executor
    private WrRestApi client
    private final Path wrapperFile
    private final Path outputFile
    private final Path errorFile
    private final Path logFile
    private final Path scriptFile
    private final Path inputFile
    private final Path traceFile
    private String jobId

    /**
     * Batch context shared between multiple task handlers
     */
    private BatchContext<String,Map> context

    WrTaskHandler(TaskRun task, WrExecutor executor) {
        super(task)
        this.executor = executor
        this.client = executor.getClient()

        this.logFile = task.workDir?.resolve(TaskRun.CMD_LOG)
        this.scriptFile = task.workDir?.resolve(TaskRun.CMD_SCRIPT)
        this.inputFile =  task.workDir?.resolve(TaskRun.CMD_INFILE)
        this.outputFile = task.workDir?.resolve(TaskRun.CMD_OUTFILE)
        this.errorFile = task.workDir?.resolve(TaskRun.CMD_ERRFILE)
        this.wrapperFile = task.workDir?.resolve(TaskRun.CMD_RUN)
        this.traceFile = task.workDir?.resolve(TaskRun.CMD_TRACE)
    }

    void batch( BatchContext<String,Map> context ) {
        if( jobId ) {
            context.collect(jobId)
            this.context = context
        }
    }

    private String jobIdsToString(Collection batchIds) {
        final MAX=10
        final sz = batchIds.size()
        batchIds.size()<=MAX ? batchIds.join(',').toString() : batchIds.take(MAX).join(',').toString() + ", ... other ${sz-MAX} omitted"
    }

    /**
     * Retrieve Batch job status information
     *
     * @param jobId The Batch job ID
     * @return The associated job details in a Map or {@code null} if no information is found
     */
    protected Map getJob(String jobId) {
        Collection batchIds
        if( context ) {
            // check if this response is cached in the batch collector
            if( context.contains(jobId) ) {
                log.trace "[wr] hit cache when getting job=$jobId"
                return context.get(jobId)
            }
            log.trace "[wr] missed cache when getting job=$jobId"
            batchIds = context.getBatchFor(jobId, 1000)
        }
        else {
            batchIds = [jobId]
        }

        // retrieve the status for the specified batch of jobs
        final String ids = jobIdsToString(batchIds)
        log.trace "[wr] getting jobs=${jobIdsToString(batchIds)}"
        List<Map> jobs = client.status(batchIds.join(',').toString())
        if( !jobs ) {
            log.debug "[wr] cannot retrieve status for job=$jobId"
            return null
        }

        Map result=null
        jobs.each {
            String id = it."Key" as String
            // cache the response in the batch collector
            context?.put( id, it )
            // return the job detail for the specified job
            if( id == jobId )
                result = it
        }
        if( !result ) {
            log.debug "[wr] cannot find status for job=$jobId"
        }

        return result
    }

    @Override
    boolean checkIfRunning() {
        if( !jobId || !isSubmitted() )
            return false
        final job = getJob(jobId)
        final result = job?.State in STARTED_STATUSES
        if( result ) {
            log.trace "[wr] Task started > $task.name"
            this.status = TaskStatus.RUNNING
        }
        return result
    }

    @Override
    boolean checkIfCompleted() {
        if( !isRunning() )
            return false
        final job = getJob(jobId)
        final done = job?.State in COMPLETE_STATUSES
        if( done ) {
            // finalize the task
            log.trace "[wr] Task completed > $task.name"
            if (job.Exited) {
                task.exitStatus = job.Exitcode as Integer
            } else {
                task.exitStatus = Integer.MAX_VALUE
            }
            task.stdout = outputFile
            task.stderr = errorFile
            status = TaskStatus.COMPLETED
        }
        return done
    }

    @Override
    void kill() {
        if( jobId )
            client.cancel(jobId)
        else
            log.trace "[wr] Invalid kill request -- missing jobId"
    }

    @Override
    void submit() {
        // this is not actually called by anything since WrMonitor does batched
        // submits direcly to the client, but this is implemented anyway
        List<Map> jobs = client.add([submitArgs()])
        jobId = jobs[0]."Key" as String
        // log.debug("[wr] submitted job $jobId")
        status = TaskStatus.SUBMITTED
    }

    List submitArgs() {
        // create task wrapper
        final bash = new WrBashBuilder(task)
        bash.build()

        WrFileCopyStrategy copyStrategy = bash.copyStrategy as WrFileCopyStrategy
        String wrapperPath = copyStrategy.wrWorkPath(wrapperFile)

        return ["/bin/bash $wrapperPath", task, copyStrategy]
    }

    void submitted(String id) {
        jobId = id
        status = TaskStatus.SUBMITTED
    }

}