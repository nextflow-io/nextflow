/*
 * Copyright 2018, WuxiNextcode
 * Copyright 2018, Centre for Genomic Regulation (CRG)
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

package nextflow.cloud.google.pipelines

import java.nio.file.Path

import com.google.api.services.genomics.v2alpha1.model.Event
import com.google.api.services.genomics.v2alpha1.model.Metadata
import com.google.api.services.genomics.v2alpha1.model.Mount
import com.google.api.services.genomics.v2alpha1.model.Operation
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.exception.ProcessUnrecoverableException
import nextflow.processor.TaskBean
import nextflow.processor.TaskHandler
import nextflow.processor.TaskRun
import nextflow.processor.TaskStatus
import nextflow.util.Escape

/**
 * Task handler for Google Pipelines.
 *
 * @author Ólafur Haukur Flygenring <olafurh@wuxinextcode.com>
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class GooglePipelinesTaskHandler extends TaskHandler {

    static private List<String> UNSTAGE_CONTROL_FILES = [TaskRun.CMD_ERRFILE, TaskRun.CMD_OUTFILE, TaskRun.CMD_LOG, TaskRun.CMD_EXIT ]

    static private String DEFAULT_INSTANCE_TYPE = 'n1-standard-1'

    GooglePipelinesExecutor executor
    
    @PackageScope GooglePipelinesConfiguration pipelineConfiguration

    private TaskBean taskBean

    private Path exitFile
    private Path wrapperFile
    private Path outputFile
    private Path errorFile
    private Path logFile
    private Path scriptFile
    private Path inputFile
    private Path stubFile
    private Path traceFile

    private String instanceType
    private String mountPath

    final static String diskName = "nf-pipeline-work"
    final static String fileCopyImage = "google/cloud-sdk:alpine"

    private Mount sharedMount

    private Operation operation

    private Metadata metadata

    @PackageScope List<String> stagingCommands = []

    @PackageScope List<String> unstagingCommands = []

    GooglePipelinesTaskHandler(TaskRun task, GooglePipelinesExecutor executor, GooglePipelinesConfiguration pipelineConfiguration) {
        super(task)

        this.executor = executor
        this.taskBean = new TaskBean(task)
        this.pipelineConfiguration = pipelineConfiguration

        this.logFile = task.workDir.resolve(TaskRun.CMD_LOG)
        this.scriptFile = task.workDir.resolve(TaskRun.CMD_SCRIPT)
        this.inputFile = task.workDir.resolve(TaskRun.CMD_INFILE)
        this.outputFile = task.workDir.resolve(TaskRun.CMD_OUTFILE)
        this.errorFile = task.workDir.resolve(TaskRun.CMD_ERRFILE)
        this.exitFile = task.workDir.resolve(TaskRun.CMD_EXIT)
        this.wrapperFile = task.workDir.resolve(TaskRun.CMD_RUN)
        this.stubFile = task.workDir.resolve(TaskRun.CMD_STUB)
        this.traceFile = task.workDir.resolve(TaskRun.CMD_TRACE)

        //Set the mount path to be the workdir that is parent of the hashed directories.
        this.mountPath = task.workDir.parent.parent.toString()

        //Get the instanceType to use for this task
        instanceType = executor.getSession().config.navigate("cloud.instanceType") ?: DEFAULT_INSTANCE_TYPE

        validateConfiguration()

        log.debug "[GPAPI] Created handler for task '${task.name}'"
    }


    /* ONLY FOR TESTING PURPOSE */
    protected GooglePipelinesTaskHandler() {

    }

    void validateConfiguration() {
        if (!task.container) {
            throw new ProcessUnrecoverableException("No container image specified for process $task.name -- Either specify the container to use in the process definition or with 'process.container' value in your config")
        }
        if (!instanceType) {
            throw new ProcessUnrecoverableException("No instance type specified for process $task.name -- Please provide a 'cloud.instanceType' definition in your config")
        }
    }

    /**
     * @return
     *      It should return {@code true} only the very first time the
     *      the task status transitions from SUBMITTED to RUNNING, in all other
     *      cases if should return {@code false}
     */
    @Override
    boolean checkIfRunning() {
        if( operation==null || !isSubmitted() )
            return false

        // note, according to the semantic of this method
        // the handler status has to be changed to RUNNING either
        // if the operation is still running or it has completed
        def result = executor.helper.checkOperationStatus(operation)
        if( result!=null ) {
            status = TaskStatus.RUNNING
            operation = result
            return true
        }
        return false
    }

    @Override
    boolean checkIfCompleted() {
        if( !isRunning() )
            return false
        
        operation = executor.helper.checkOperationStatus(operation)

        def events = extractRuntimeDataFromOperation()
        events?.reverse()?.each {
            log.trace "[GPAPI] New event for task '$task.name' - time: ${it.get("timestamp")} - ${it.get("description")}"
        }

        if (operation.getDone()) {
            log.debug "[GPAPI] Task '$task.name' complete. Start Time: ${metadata?.getStartTime()} - End Time: ${metadata?.getEndTime()}"

            // finalize the task
            Integer xs = readExitFile()
            //Use the status from the exitStatus file if it exists. Else use the exit code from the pipeline operation
            task.stdout = outputFile
            task.exitStatus = xs != null ? xs :  operation.getError()?.getCode()
            task.stderr = xs  != null ?  errorFile : operation.getError()?.getMessage()

            /* If we get error 10 or 14 from the pipelines api we'll retry
             * since it means that the vm instance was probably preemptied.
             *
             * This retry method is a bit crude, find a cleaner method if it is available
             */
            if(!xs && (task.exitStatus == 10 || task.exitStatus == 14)) {
                task.config.setProperty("maxRetries",task.failCount+1)
                task.config.setProperty("errorStrategy","retry")
            }
            status = TaskStatus.COMPLETED
            return true
        } else
            return false
    }

    private List<Event> extractRuntimeDataFromOperation() {
        def metadata = operation.getMetadata() as Metadata
        if (!this.metadata) {
            this.metadata = metadata
            return metadata?.getEvents()
        } else {
            //Get the new events
            def delta = metadata.getEvents().size() - this.metadata.getEvents().size()
            this.metadata = metadata
            return delta > 0 ? metadata.getEvents().take(delta) : [] as List<Event>
        }
    }

    @PackageScope Integer readExitFile() {
        try {
            exitFile.text as Integer
        }
        catch (Exception e) {
            log.debug "[GPAPI] Cannot read exitstatus for task: `$task.name`", e
            null
        }
    }

    @Override
    void kill() {
        if( !operation ) return
        log.debug "[GPAPI] Killing pipeline '${operation.name}'"
        executor.helper.cancelOperation(operation)
    }

    @Override
    void submit() {
        createTaskWrapper()
        final req = createPipelineRequest()
        operation = submitPipeline(req)
        status = TaskStatus.SUBMITTED
        log.trace "[GPAPI] Submitted task '$task.name. Assigned Pipeline operation name = '${operation.getName()}'"
    }

    @PackageScope
    void createTaskWrapper() {
        new GooglePipelinesScriptLauncher(this.taskBean, this) .build()
    }

    @PackageScope
    Operation submitPipeline(GooglePipelinesSubmitRequest request) {
        executor.helper.submitPipeline(request)
    }

    @PackageScope
    GooglePipelinesSubmitRequest createPipelineRequest() {
        String stagingScript = "mkdir -p ${task.workDir.toString()}"
        if( stagingCommands )
            stagingScript += "; (${stagingCommands.join("; ")}) 2>&1 > $task.workDir/${TaskRun.CMD_LOG}"

        String mainScript = "cd $task.workDir; bash ${TaskRun.CMD_RUN} 2>&1 | tee -a ${TaskRun.CMD_LOG}"

        /*
         * -m = run in parallel
         * -q = quiet mode
         * cp = copy
         * -P = preserve POSIX attributes
         * -c = continues on errors
         * -r = recursive copy
         */
        def gsCopyPrefix = "gsutil -m -q cp -c -P"

        //Copy the logs provided by Google Pipelines for the pipeline to our work dir.
        def unstaging = []
        unstaging << "[[ \$GOOGLE_PIPELINE_FAILED == 1 ]] && $gsCopyPrefix -r /google/ ${task.workDir.toUriString()} || true".toString()
        unstaging.addAll(unstagingCommands)

        //add the task output files to unstaging command list
        for( String it : UNSTAGE_CONTROL_FILES ) {
            unstaging << "$gsCopyPrefix $task.workDir/${Escape.path(it)} ${task.workDir.toUriString()} || true".toString()
        }

        //Create the mount for out work files.
        sharedMount = executor.helper.configureMount(diskName, mountPath)

        def req = new GooglePipelinesSubmitRequest()
        req.instanceType = instanceType
        req.project = pipelineConfiguration.project
        req.zone = pipelineConfiguration.zone
        req.region = pipelineConfiguration.region
        req.diskName = diskName
        req.preemptible = pipelineConfiguration.preemptible
        req.taskName = "nf-$task.hash"
        req.containerImage = task.container
        req.fileCopyImage = fileCopyImage
        req.stagingScript = stagingScript
        req.mainScript = mainScript
        req.unstagingScript = unstaging.join("; ").leftTrim()
        req.sharedMount = sharedMount
        return req
    }


}