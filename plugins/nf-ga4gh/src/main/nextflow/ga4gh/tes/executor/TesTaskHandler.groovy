/*
 * Copyright 2013-2024, Seqera Labs
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

package nextflow.ga4gh.tes.executor

import nextflow.ga4gh.tes.client.model.TesFileType
import nextflow.processor.TaskBean

import static nextflow.processor.TaskStatus.COMPLETED
import static nextflow.processor.TaskStatus.RUNNING

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.exception.ProcessUnrecoverableException
import nextflow.executor.BashWrapperBuilder
import nextflow.ga4gh.tes.client.api.TaskServiceApi
import nextflow.ga4gh.tes.client.model.TesExecutor as TesExecutorModel
import nextflow.ga4gh.tes.client.model.TesInput
import nextflow.ga4gh.tes.client.model.TesOutput
import nextflow.ga4gh.tes.client.model.TesResources
import nextflow.ga4gh.tes.client.model.TesState
import nextflow.ga4gh.tes.client.model.TesTask
import nextflow.processor.TaskConfig
import nextflow.processor.TaskHandler
import nextflow.processor.TaskRun
import nextflow.processor.TaskStatus
import nextflow.script.params.FileOutParam
import nextflow.util.Escape
import nextflow.util.MemoryUnit
/**
 * Handle execution phases for a task executed by a TES executor
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class TesTaskHandler extends TaskHandler {

    static final private String WORK_DIR = '/work'

    final List<TesState> COMPLETE_STATUSES = [TesState.COMPLETE, TesState.EXECUTOR_ERROR, TesState.SYSTEM_ERROR, TesState.CANCELED]

    final List<TesState> STARTED_STATUSES = [TesState.INITIALIZING, TesState.RUNNING, TesState.PAUSED] + COMPLETE_STATUSES

    private TesExecutor executor

    private final Path exitFile

    private final Path wrapperFile

    private final Path outputFile

    private final Path errorFile

    private final Path logFile

    private final Path scriptFile

    private final Path inputFile

    private final Path traceFile

    private TaskServiceApi client

    private String requestId

    /** only for testing purpose -- do not use */
    protected TesTaskHandler() {}

    TesTaskHandler(TaskRun task, TesExecutor executor) {
        super(task)
        this.executor = executor
        this.client = executor.getClient()

        this.logFile = task.workDir.resolve(TaskRun.CMD_LOG)
        this.scriptFile = task.workDir.resolve(TaskRun.CMD_SCRIPT)
        this.inputFile =  task.workDir.resolve(TaskRun.CMD_INFILE)
        this.outputFile = task.workDir.resolve(TaskRun.CMD_OUTFILE)
        this.errorFile = task.workDir.resolve(TaskRun.CMD_ERRFILE)
        this.exitFile = task.workDir.resolve(TaskRun.CMD_EXIT)
        this.wrapperFile = task.workDir.resolve(TaskRun.CMD_RUN)
        this.traceFile = task.workDir.resolve(TaskRun.CMD_TRACE)
    }

    protected String getRequestId() { requestId }

    @Override
    boolean checkIfRunning() {

        if( requestId && isSubmitted() ) {
            final response = client.getTask(requestId, null)
            final started = response.state in STARTED_STATUSES
            if( started ) {
                log.trace "[TES] Task started > $task.name"
                status = RUNNING
                return true
            }
        }

        return false
    }

    @Override
    boolean checkIfCompleted() {
        if( !isRunning() ) {
            return false
        }

        final response = client.getTask(requestId, null)
        if( response.state in COMPLETE_STATUSES ) {
            // finalize the task
            log.trace "[TES] Task completed > $task.name"
            task.exitStatus = readExitStatus()
            task.stdout = outputFile
            task.stderr = errorFile
            status = COMPLETED
            return true
        }

        return false
    }

    private int readExitStatus() {
        try {
            return client
                .getTask(requestId, 'FULL')
                .logs[0]
                .logs[0]
                .exitCode
        }
        catch( Exception e ) {
            log.trace "[TES] Cannot read exitstatus for task: `$task.name` | ${e.message}"
            return Integer.MAX_VALUE
        }
    }

    @Override
    void kill() {
        if( requestId )
            client.cancelTask(requestId)
        else
            log.trace "[TES] Invalid kill request -- missing requestId"
    }

    @Override
    void submit() {

        // create task wrapper
        String remoteBinDir = executor.getRemoteBinDir()
        final bash = newTesBashBuilder(task, remoteBinDir)
        bash.build()

        final body = newTesTask()

        // submit the task
        final task = client.createTask(body)
        requestId = task.id
        status = TaskStatus.SUBMITTED
    }

    protected TesBashBuilder newTesBashBuilder(TaskRun task, String remoteBinDir) {
        final builder = new TesBashBuilder(task, remoteBinDir)
        builder.headerScript = "NXF_CHDIR=${Escape.path(WORK_DIR)}"
        return builder
    }

    protected TesTask newTesTask() {
        if( !task.container )
            throw new ProcessUnrecoverableException("Process `${task.lazyName()}` failed because the container image was not specified")

        final exec = new TesExecutorModel()
        exec.command = List.of('bash', '-o', 'pipefail', '-c', "bash ${WORK_DIR}/${TaskRun.CMD_RUN} 2>&1 | tee ${WORK_DIR}/${TaskRun.CMD_LOG}".toString())
        exec.image = task.container
        exec.workdir = WORK_DIR

        final body = new TesTask()

        // add task control files
        body.addInputsItem(inItem(scriptFile))
        body.addInputsItem(inItem(wrapperFile))

        for( Path path : executor.getRemoteBinFiles() )
            body.addInputsItem(inItemFromBin(path))

        // add task input files
        if( inputFile.exists() )
            body.addInputsItem(inItem(inputFile))

        task.getInputFilesMap().each { String name, Path path ->
            body.addInputsItem(inItem(path, name))
        }

        // add the task output files
        body.addOutputsItem(outItem(outputFile.name))
        body.addOutputsItem(outItem(errorFile.name))
        body.addOutputsItem(outItem(logFile.name))
        body.addOutputsItem(outItem(exitFile.name))

        // set requested resources
        body.setResources(getResources(task.config))

        for( final param : task.getOutputsByType(FileOutParam).keySet() ) {
            final type = param.type == 'dir'
                ? 'DIRECTORY'
                : 'FILE'
            for( final pattern : param.getFilePatterns(task.context, task.workDir) )
                body.addOutputsItem(outItem(pattern, type))
        }

        body.setName(task.getName())

        // add the executor
        body.executors = List.of(exec)

        return body
    }

    private TesResources getResources(TaskConfig config) {
        final res = new TesResources()
            .cpuCores(config.getCpus())
            .ramGb(toGiga(config.getMemory()))
            .diskGb(config.getDisk()?.toGiga())
        log.trace("[TES] Adding resource request: $res")
        // @TODO preemptible
        // @TODO zones
        return res
    }

    private Double toGiga(MemoryUnit size) {
        // 1073741824 = 1GB
        return size != null ? ((double)size.bytes)/1073741824 : null
    }

    private TesInput inItem(Path realPath, String fileName = null) {
        final result = new TesInput()
        result.url = realPath.toUriString()
        result.path = fileName ? "$WORK_DIR/$fileName" : "$WORK_DIR/${realPath.getName()}"
        result.type = realPath.isDirectory() ? 'DIRECTORY' : 'FILE'

        // fix local path for TES Azure
        result.url = fixTesAzureLocalPath(result.url)

        log.trace "[TES] Adding INPUT file: $result"
        return result
    }

    private TesInput inItemFromBin(Path realPath) {
        final result = new TesInput()
        result.url = realPath.toUriString()
        result.path = realPath.toString()
        result.type = 'FILE'

        // fix local path for TES Azure
        result.url = fixTesAzureLocalPath(result.url)

        log.trace "[TES] Adding INPUT file: $result"
        return result
    }

    private TesOutput outItem( String fileName, String type = null ) {
        final result = new TesOutput()
        if( fileName.contains('*') || fileName.contains('?') ) {
            result.path = "$WORK_DIR/$fileName"
            result.pathPrefix = WORK_DIR
            result.url = task.workDir.toUriString()
        }
        else {
            result.path = "$WORK_DIR/$fileName"
            result.url = task.workDir.resolve(fileName).toUriString()
        }
        if( type != null )
            result.type = type

        // fix local path for TES Azure
        result.url = fixTesAzureLocalPath(result.url)

        log.trace "[TES] Adding OUTPUT file: $result"
        return result
    }

    /**
     * Fix local paths when using TES Azure, which requires the
     * following AZ path:
     *
     *   az://<blob>/<path>
     *
     * to be formatted as:
     *
     *   /<storage-account>/<blob>/<path>
     *
     * @param url
     */
    private String fixTesAzureLocalPath(String url) {
        if( !url.startsWith('az://') )
            return url
        final storageAccount = executor.getAzureStorageAccount()
        if( !executor.getEndpoint().contains('azure.com') || !storageAccount )
            return url
        return url.replaceAll('az://', "/${storageAccount}/")
    }

}
