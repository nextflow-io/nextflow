/*
 * Copyright 2018, WuxiNextcode
 * Copyright 2018-2019, Centre for Genomic Regulation (CRG)
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
import groovy.json.JsonOutput
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.exception.ProcessSubmitException
import nextflow.exception.ProcessUnrecoverableException
import nextflow.processor.TaskBean
import nextflow.processor.TaskHandler
import nextflow.processor.TaskRun
import nextflow.processor.TaskStatus
import nextflow.trace.TraceRecord
import nextflow.util.Escape
import nextflow.util.MemoryUnit
/**
 * Task handler for Google Pipelines.
 *
 * @author Ólafur Haukur Flygenring <olafurh@wuxinextcode.com>
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class GooglePipelinesTaskHandler extends TaskHandler {

    static private List<String> UNSTAGE_CONTROL_FILES = [TaskRun.CMD_ERRFILE, TaskRun.CMD_OUTFILE, TaskRun.CMD_LOG, TaskRun.CMD_EXIT ]

    GooglePipelinesExecutor executor
    
    @PackageScope GooglePipelinesConfiguration pipelineConfiguration

    private TaskBean taskBean

    private Path exitFile

    private Path outputFile

    private Path errorFile

    private String machineType

    private String mountPath

    final static String diskName = "nf-pipeline-work"

    final static String fileCopyImage = "google/cloud-sdk:alpine"

    private Mount sharedMount

    private Operation operation

    private Metadata metadata

    private String pipelineId

    @PackageScope List<String> stagingCommands = []

    @PackageScope List<String> unstagingCommands = []

    GooglePipelinesTaskHandler(TaskRun task, GooglePipelinesExecutor executor, GooglePipelinesConfiguration pipelineConfiguration) {
        super(task)

        this.executor = executor
        this.taskBean = new TaskBean(task)
        this.pipelineConfiguration = pipelineConfiguration

        this.outputFile = task.workDir.resolve(TaskRun.CMD_OUTFILE)
        this.errorFile = task.workDir.resolve(TaskRun.CMD_ERRFILE)
        this.exitFile = task.workDir.resolve(TaskRun.CMD_EXIT)

        //Set the mount path to be the workdir that is parent of the hashed directories.
        this.mountPath = task.workDir.parent.parent.toString()

        //Get the machineType to use for this task
        this.machineType = getProcessMachineType()

        validateConfiguration()
    }


    /* ONLY FOR TESTING PURPOSE */
    protected GooglePipelinesTaskHandler() {

    }

    void validateConfiguration() {
        if (!task.container) {
            throw new ProcessUnrecoverableException("No container image specified for process $task.name -- Either specify the container to use in the process definition or with 'process.container' value in your config")
        }
    }

    //
    // Use the process cpus and memory directives to create custom GCP instance
    // If memory not specified, default to 1 GB per cpu.  An absence of the cpus directive defaults to 1 cpu.
    // If the process machineType is defined, use that instead of cpus/memory (must be a predefined GCP machine type)
    // The deprecated cloud.instanceType will be used if specified and machineType is not specified
    String getProcessMachineType() {
        String machineType = getProcessMachineType(executor.getSession().config.navigate("cloud.instanceType"), task.config.machineType, task.config.cpus, task.config.memory)

        log.trace "[GPAPI] Task: $task.name - Instance Type: $machineType"

        return machineType
    }

    String getProcessMachineType(String cloudInstanceType, String taskMachineType, int cpus, MemoryUnit memory) {

        if (cloudInstanceType != null) {
            log.warn1("Configuration setting [cloud.instanceType] is deprecated.  Please use either the process cpus/memory directives or process.machineType instead.", firstOnly: true, cacheKey: GooglePipelinesTaskHandler.class)
        }

        String machineType = taskMachineType ?: cloudInstanceType

        if (machineType == null) {
            long megabytes = memory != null ? memory.mega : cpus*1024
            machineType = 'custom-' + cpus + '-' + megabytes
        }

        return machineType
    }


    protected void logEvents(Operation operation) {
        final events = getEventsFromOp(operation)
        if( !events )
            return

        final warns = new HashSet()
        for( Event e : events ) {
            def d = e.getDescription()
            if( d?.contains('resource_exhausted') )
                warns << d
        }

        if( warns ) {
            log.debug "[GPAPI] New event > $task.name - Pipeline Id: $pipelineId\n${prettyPrint(events)}"
            for( String w : warns ) log.warn1("Google Pipelines > $w")
        }
        else if( log.isTraceEnabled() ) {
            log.trace "[GPAPI] New event > $task.name - Pipeline Id: $pipelineId\n${prettyPrint(events)}"
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
        final result = executor.helper.checkOperationStatus(operation)
        if( result == null )
            return false
        logEvents(result)

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

        final resultOp = executor.helper.checkOperationStatus(operation)
        if( !resultOp )
            return false

        operation = resultOp
        logEvents(operation)

        if (operation.getDone()) {
            log.debug "[GPAPI] Task complete > $task.name - Start Time: ${metadata?.getStartTime()} - End Time: ${metadata?.getEndTime()}"

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

    @PackageScope
    List<Event> getEventsFromOp(Operation operation) {
        final metadata = (Metadata)operation.getMetadata()
        List<Event> result
        if (!this.metadata) {
            this.metadata = metadata
            result = metadata != null ? metadata.getEvents() : Collections.<Event>emptyList()
        }
        else {
            //Get the new events
            def delta = metadata.getEvents().size() - this.metadata.getEvents().size()
            this.metadata = metadata
            result = delta > 0 ? metadata.getEvents().take(delta) : Collections.<Event>emptyList()
        }
        return result.reverse()
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
        log.debug "[GPAPI] Killing task > $task.name - Pipeline Id: $pipelineId"
        executor.helper.cancelOperation(operation)
    }

    @Override
    void submit() {
        createTaskWrapper()
        final req = createPipelineRequest()
        log.trace "[GPAPI] Task created > $task.name - Request: $req"

        operation = submitPipeline(req)
        if( operation == null )
            throw new ProcessSubmitException("Failed to submit task with name: $task.name")
        pipelineId = getPipelineIdFromOp(operation)
        status = TaskStatus.SUBMITTED

        if( log.isTraceEnabled() ) {
            log.trace "[GPAPI] Task submitted > $task.name - Pipeline Id: $pipelineId; Operation:\n${prettyPrint(operation)}"
        }
        else {
            log.debug "[GPAPI] Task submitted > $task.name - Pipeline Id: $pipelineId"
        }
    }

    @PackageScope
    String getPipelineIdFromOp(Operation operation) {
        assert operation?.getName()
        operation.getName().tokenize('/')[-1]
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

        // Copy the logs provided by Google Pipelines for the pipeline to our work dir.
        def unstaging = []
        unstaging << "cd $task.workDir"
        unstaging << '[[ $GOOGLE_PIPELINE_FAILED == 1 || $NXF_DEBUG ]] && ' + unstage('/google/')
        unstaging.addAll(unstagingCommands)
        // unstage trace file
        unstaging << "[[ -f $TaskRun.CMD_TRACE ]] && ${unstage(TaskRun.CMD_TRACE)}"
        // add the task output files to unstaging command list
        for( String it : UNSTAGE_CONTROL_FILES ) {
            unstaging << unstage(it)
        }

        //Create the mount for out work files.
        sharedMount = executor.helper.configureMount(diskName, mountPath)

        def req = new GooglePipelinesSubmitRequest()
        req.machineType = machineType
        req.project = pipelineConfiguration.project
        req.zone = pipelineConfiguration.zone
        req.region = pipelineConfiguration.region
        req.diskName = diskName
        req.diskSizeGb = task.config.disk?.giga
        req.preemptible = pipelineConfiguration.preemptible
        req.taskName = "nf-$task.hash"
        req.containerImage = task.container
        req.fileCopyImage = fileCopyImage
        req.stagingScript = stagingScript
        req.mainScript = mainScript
        req.unstagingScript = unstaging.join("; ").trim()
        req.sharedMount = sharedMount
        return req
    }

    @PackageScope
    String unstage(String local){
        /*
         * -m = run in parallel
         * -q = quiet mode
         * cp = copy
         * -R = recursive copy
         */
        "gsutil -m -q cp -R ${Escape.path(local)} ${task.workDir.toUriString()} || true"
    }

    TraceRecord getTraceRecord() {
        def result = super.getTraceRecord()
        result.put('native_id', pipelineId)
        return result
    }

    String prettyPrint(Operation op) {
        JsonOutput.prettyPrint( JsonOutput.toJson(op) )
    }

    static String prettyPrint(List<Event> events) {
        JsonOutput.prettyPrint( JsonOutput.toJson(events) )
    }

}
