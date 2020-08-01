/*
 * Copyright 2019, Google Inc
 * Copyright 2018, WuxiNextcode
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

package nextflow.cloud.google.lifesciences

import static nextflow.cloud.google.lifesciences.GoogleLifeSciencesHelper.*

import java.nio.file.Path

import com.google.api.services.lifesciences.v2beta.model.Event
import com.google.api.services.lifesciences.v2beta.model.Metadata
import com.google.api.services.lifesciences.v2beta.model.Mount
import com.google.api.services.lifesciences.v2beta.model.Operation
import groovy.json.JsonOutput
import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.cloud.types.CloudMachineInfo
import nextflow.cloud.types.PriceModel
import nextflow.exception.ProcessSubmitException
import nextflow.exception.ProcessUnrecoverableException
import nextflow.processor.TaskBean
import nextflow.processor.TaskHandler
import nextflow.processor.TaskRun
import nextflow.processor.TaskStatus
import nextflow.trace.TraceRecord
import nextflow.util.MemoryUnit
/**
 * Task handler for Google Pipelines.
 *
 * @author Ólafur Haukur Flygenring <olafurh@wuxinextcode.com>
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@Slf4j
class GoogleLifeSciencesTaskHandler extends TaskHandler {

    public final static String DEFAULT_DISK_NAME = "nf-pipeline-work"

    GoogleLifeSciencesExecutor executor

    private TaskBean taskBean

    private Path exitFile

    private Path outputFile

    private Path errorFile

    private Operation operation

    private Metadata metadata

    private String pipelineId

    private GoogleLifeSciencesHelper helper

    private volatile String assignedZone

    private volatile String assignedInstance

    GoogleLifeSciencesTaskHandler(TaskRun task, GoogleLifeSciencesExecutor executor) {
        super(task)

        this.executor = executor
        this.taskBean = new TaskBean(task)
        this.helper = executor.helper

        this.outputFile = task.workDir.resolve(TaskRun.CMD_OUTFILE)
        this.errorFile = task.workDir.resolve(TaskRun.CMD_ERRFILE)
        this.exitFile = task.workDir.resolve(TaskRun.CMD_EXIT)

        validateConfiguration()
    }


    /* ONLY FOR TESTING PURPOSE */
    protected GoogleLifeSciencesTaskHandler() {

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
    protected String getMachineType() {
        String machineType = getMachineType0(task.config.getMachineType(), task.config.getCpus(), task.config.getMemory())
        log.trace "[GLS] Task: $task.name - Instance Type: $machineType"
        return machineType
    }

    protected String getMachineType0(String taskMachineType, int cpus, MemoryUnit memory) {

        String machineType = taskMachineType

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
        for( Event ev : events ) {
            log.trace "[GLS] task $task.name > event=$ev"
            final d = ev.getDescription()
            // collect warnings
            if( d?.contains('resource_exhausted') )
                warns << d
            // fetch the assigned zone
            if( !assignedZone && ev.getWorkerAssigned() )
                assignedZone = ev.getWorkerAssigned().getZone()
            // fetch the assigned instance
            if( !assignedInstance && ev.getWorkerAssigned() )
                assignedInstance = ev.getWorkerAssigned().getInstance()
            // dump SSH daemon info
            if( ev.getContainerStarted() && d?.contains(SSH_DAEMON_NAME) )
                log.debug "[GLS] SSH daemon IP ${ev.getContainerStarted().getIpAddress()}; connect command: `gcloud compute --project ${executor.config.project} ssh --zone ${assignedZone} ${assignedInstance}`"
        }

        if( warns ) {
            log.debug "[GLS] New event > $task.name - Pipeline Id: $pipelineId\n${prettyPrint(events)}"
            for( String w : warns ) log.warn1("Google Pipelines > $w")
        }
        else if( log.isTraceEnabled() ) {
            log.trace "[GLS] New event > $task.name - Pipeline Id: $pipelineId\n${prettyPrint(events)}"
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
        final result = helper.checkOperationStatus(operation)
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

        final resultOp = helper.checkOperationStatus(operation)
        if( !resultOp )
            return false

        operation = resultOp
        logEvents(operation)

        if (operation.getDone()) {
            log.debug "[GLS] Task complete > $task.name - Start Time: ${metadata?.getStartTime()} - End Time: ${metadata?.getEndTime()}"

            // finalize the task
            Integer xs = readExitFile()
            //Use the status from the exitStatus file if it exists. Else use the exit code from the pipeline operation
            task.stdout = outputFile
            task.exitStatus = xs != null ? xs :  operation.getError()?.getCode()
            task.stderr = xs  != null ?  errorFile : operation.getError()?.getMessage()
            status = TaskStatus.COMPLETED
            return true
        } else
            return false
    }

    @PackageScope
    List<Event> getEventsFromOp(Operation operation) {
        final metadata = (Metadata)operation.getMetadata()
        if( !metadata?.getEvents() )
            return Collections.<Event>emptyList()

        List<Event> result
        if ( !this.metadata ) {
            this.metadata = metadata
            result = metadata.getEvents()
        }
        else {
            //Get the new events
            def delta = metadata.getEvents().size() - this.metadata.getEvents().size()
            this.metadata = metadata
            result = delta > 0 ? metadata.getEvents().take(delta) : null
        }
        return result != null ? result.reverse() : Collections.<Event>emptyList()
    }

    @PackageScope Integer readExitFile() {
        try {
            exitFile.text as Integer
        }
        catch (Exception e) {
            log.debug "[GLS] Cannot read exitstatus for task: `$task.name` | ${e.message}"
            null
        }
    }

    @Override
    void kill() {
        if( !operation ) return
        log.debug "[GLS] Killing task > $task.name - Pipeline Id: $pipelineId"
        helper.cancelOperation(operation)
    }

    @Override
    void submit() {
        createTaskWrapper()
        final req = createPipelineRequest()
        log.trace "[GLS] Task created > $task.name - Request: $req"

        operation = submitPipeline(req)
        if( operation == null )
            throw new ProcessSubmitException("Failed to submit task with name: $task.name")
        pipelineId = getPipelineIdFromOp(operation)
        status = TaskStatus.SUBMITTED

        if( log.isTraceEnabled() ) {
            log.trace "[GLS] Task submitted > $task.name - Pipeline Id: $pipelineId; Operation:\n${prettyPrint(operation)}"
        }
        else {
            log.debug "[GLS] Task submitted > $task.name - Pipeline Id: $pipelineId"
        }
    }

    @PackageScope
    String getPipelineIdFromOp(Operation operation) {
        assert operation?.getName()
        operation.getName().tokenize('/')[-1]
    }

    @PackageScope
    void createTaskWrapper() {
        new GoogleLifeSciencesScriptLauncher(this.taskBean, this) .build()
    }

    @PackageScope
    Operation submitPipeline(GoogleLifeSciencesSubmitRequest request) {
        helper.submitPipeline(request)
    }

    @PackageScope
    GoogleLifeSciencesSubmitRequest createPipelineRequest() {
        //Create the mount for out work files.
        def req = new GoogleLifeSciencesSubmitRequest()
        req.machineType = getMachineType()
        req.project = executor.config.project
        req.zone = executor.config.zones
        req.region = executor.config.regions
        req.diskName = DEFAULT_DISK_NAME
        req.diskSizeGb = task.config.getDisk()?.getGiga() as Integer
        req.preemptible = executor.config.preemptible
        req.taskName = "nf-$task.hash"
        req.containerImage = task.container
        req.workDir = task.workDir
        req.sharedMount = configureMount(DEFAULT_DISK_NAME, task.workDir.toString())
        req.accelerator = task.config.getAccelerator()
        req.location = executor.config.location
        req.cpuPlatform = executor.config.cpuPlatform
        req.bootDiskSizeGb = executor.config.bootDiskSize?.toGiga() as Integer
        req.entryPoint = task.config.getContainerOptionsMap().getOrDefault('entrypoint', GoogleLifeSciencesConfig.DEFAULT_ENTRY_POINT)
        req.usePrivateAddress = executor.config.usePrivateAddress
        return req
    }

    protected Mount configureMount(String diskName, String mountPath, boolean readOnly = false) {
        new Mount().setDisk(diskName).setPath(mountPath).setReadOnly(readOnly)
    }

    TraceRecord getTraceRecord() {
        def result = super.getTraceRecord()
        result.put('native_id', pipelineId)
        result.machineInfo = getMachineInfo()
        return result
    }

    private CloudMachineInfo getMachineInfo() {
        final price = executor.config.preemptible ? PriceModel.spot : PriceModel.standard
        final result = new CloudMachineInfo(machineType, assignedZone, price)
        log.trace "[GLS] Task: $task.name > cloud-info=$result"
        return result
    }

    static String prettyPrint(Operation op) {
        JsonOutput.prettyPrint( JsonOutput.toJson(op) )
    }

    static String prettyPrint(List<Event> events) {
        JsonOutput.prettyPrint( JsonOutput.toJson(events) )
    }

}
