package nextflow.prov.openlineage


import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import io.openlineage.client.OpenLineage
import nextflow.Session
import nextflow.processor.TaskRun
import nextflow.trace.TraceRecord

import java.nio.file.Path
import java.nio.file.Paths
import java.time.ZonedDateTime
import java.util.concurrent.atomic.AtomicInteger

/**
 * Factory for OpenLineage events
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Slf4j
@CompileStatic
class OpenLineageEventFactory {
    private static final OpenLineage OL = new OpenLineage(URI.create("nextflow"))
    private static final AtomicInteger publishCount = new AtomicInteger()
    private OpenLineage.Job workflow
    private OpenLineage.Run workflowRun

    OpenLineageEventFactory(){ }

    public void setWorkflowJobFromSession(Session session) {
        final namespace = session.workflowMetadata.repository ? session.workflowMetadata.repository : session.workflowMetadata.scriptFile.parent.toUri().toString()
        def name = session.workflowMetadata.manifest ? session.workflowMetadata.manifest.name : session.workflowMetadata.scriptFile.name
        name = name ? name : "main.nf"
        final workflowType = OL.newJobTypeJobFacet("BATCH", "Nextflow", "DAG")
        this.workflow = OL.newJobBuilder().namespace(namespace).name(name)
                .facets(OL.newJobFacetsBuilder().jobType(workflowType).build())
                .build()
        this.workflowRun = OL.newRunBuilder().runId(session.uniqueId).build()
        log.debug("[Openlineage] Workflow $namespace.$name job created.")
    }

    public void setWorkflowRun(UUID runId){
        this.workflowRun = OL.newRunBuilder().runId(runId).build()
    }

    public OpenLineage.RunEvent createWorkflowCompleteEvent(){
        return OL.newRunEventBuilder()
            .eventType(OpenLineage.RunEvent.EventType.COMPLETE)
            .eventTime(ZonedDateTime.now()).job(workflow).run(workflowRun)
            .build()

    }

    public OpenLineage.RunEvent createWorkflowStartEvent(){
        return OL.newRunEventBuilder().eventTime(ZonedDateTime.now())
            .eventType(OpenLineage.RunEvent.EventType.START)
            .run(this.workflowRun).job(this.workflow)
            .build()
    }

    public OpenLineage.RunEvent createWorkflowFailedEvent(){
        OL.newRunEventBuilder()
            .eventType(OpenLineage.RunEvent.EventType.FAIL)
            .eventTime(ZonedDateTime.now()).job(workflow).run(workflowRun)
            .build()
    }

    public OpenLineage.RunEvent createTaskStartEvent(TaskRun task, TraceRecord trace){
        final job = generateTaskJob(task)
        final jobRun = generateInitialTaskRun(task)
        final inputs = getInputsFromTask(task)
        return OL.newRunEventBuilder()
            .eventType(OpenLineage.RunEvent.EventType.START)
            .eventTime(ZonedDateTime.now()).job(job).run(jobRun).inputs(inputs).build()

    }

    public OpenLineage.RunEvent createTaskCompleteEvent(TaskRun task, TraceRecord trace){
        final job = generateTaskJob(task)
        final parent = generateParentFacet()
        if( task.isSuccess() ){
            return generateSuccessEvent(task, job, parent)
        } else {
            return generateFailEvent(task, job, parent)

        }
    }

    public OpenLineage.RunEvent createTaskCachedEvent(TaskRun task, TraceRecord trace){
        final job = generateTaskJob(task)
        final parent = generateParentFacet()
        return generateCachedEvent(task, job, parent)
    }

    public OpenLineage.RunEvent createPublishEvent(Path source, Path destination){
        final job = generatePublishJob()
        final parent = generateParentFacet()
        final jobRun = OL.newRun(UUID.randomUUID(),
            OL.newRunFacetsBuilder().parent(parent).build())
        return OL.newRunEventBuilder()
            .eventType(OpenLineage.RunEvent.EventType.COMPLETE)
            .eventTime(ZonedDateTime.now()).job(job).run(jobRun)
            .inputs(generateInputs([source])).outputs(generateOutputs([destination])).build()
    }


    private static OpenLineage.Job generatePublishJob(){
        final jobType = OL.newJobTypeJobFacet("BATCH", "Nextflow", "JOB")
        return OL.newJobBuilder().namespace("io.nextflow").name("publish (${publishCount.getAndIncrement()})")
            .facets(OL.newJobFacetsBuilder().jobType(jobType).build())
            .build()
    }

    private OpenLineage.Job generateTaskJob(TaskRun task){
        final jobType = OL.newJobTypeJobFacet("BATCH", "Nextflow", "JOB")
        return OL.newJobBuilder().namespace(this.workflow.namespace).name(task.name)
            .facets(OL.newJobFacetsBuilder().jobType(jobType).build())
            .build()
    }

    private OpenLineage.Run generateInitialTaskRun(TaskRun task){
        final parentRun =  generateParentFacet()
        final descFacet = generateRunFacet([hashcode: task.hash.toString(), id: task.id])
        final jobRun = OL.newRun(UUID.nameUUIDFromBytes(task.hash.asBytes()),
            OL.newRunFacetsBuilder().parent(parentRun).put("description", descFacet).build())
        return jobRun
    }

    private OpenLineage.ParentRunFacet generateParentFacet() {
        return OL.newParentRunFacet(OL.newParentRunFacetRun(this.workflowRun.getRunId()),
            OL.newParentRunFacetJob(this.workflow.namespace, this.workflow.name))
    }

    private static OpenLineage.RunEvent generateSuccessEvent(TaskRun task, OpenLineage.Job job, OpenLineage.ParentRunFacet parent){
        final jobRun = OL.newRunBuilder().runId(UUID.nameUUIDFromBytes(task.hash.asBytes()))
                .facets(OL.newRunFacetsBuilder().parent(parent).build())
                .build()
        final outputs = getOutputsFromTask(task)
        return OL.newRunEventBuilder()
                .eventType(OpenLineage.RunEvent.EventType.COMPLETE)
                .eventTime(ZonedDateTime.now())
                .job(job).run(jobRun).outputs(outputs)
                .build()
    }

    private static OpenLineage.RunEvent generateFailEvent(TaskRun task, OpenLineage.Job job, OpenLineage.ParentRunFacet parent){
        def err = task.dumpStderr()
        err = err ? err.join('\n'): ""
        final jobRun = OL.newRunBuilder().runId(UUID.nameUUIDFromBytes(task.hash.asBytes()))
                .facets(OL.newRunFacetsBuilder()
                        .errorMessage(OL.newErrorMessageRunFacet("Task $task.name failed", "Nextflow", err))
                        .parent(parent)
                        .build())
                .build()
        return OL.newRunEventBuilder()
                .eventType(OpenLineage.RunEvent.EventType.FAIL)
                .eventTime(ZonedDateTime.now()).job(job).run(jobRun)
                .build()
    }

    private static OpenLineage.RunEvent generateCachedEvent(TaskRun task, OpenLineage.Job job, OpenLineage.ParentRunFacet parent){
        final jobRun = OL.newRunBuilder().runId(UUID.nameUUIDFromBytes(task.hash.asBytes()))
                .facets(OL.newRunFacetsBuilder().parent(parent).put("abort", generateRunFacet([reason: "cached"])).build())
                .build()
        return OL.newRunEventBuilder()
                .eventType(OpenLineage.RunEvent.EventType.ABORT)
                .eventTime(ZonedDateTime.now()).job(job).run(jobRun)
                .build()
    }

    private static OpenLineage.RunFacet generateRunFacet(Map props){
        OpenLineage.RunFacet facet = OL.newRunFacet()
        props.each {facet.getAdditionalProperties().put(it.key.toString(), it.value)}
        return facet
    }

    private static List<OpenLineage.InputDataset> getInputsFromTask(TaskRun task){
        final inputs = new LinkedList<OpenLineage.InputDataset>()
        task.getInputFilesMap().each{ name,path -> inputs.add(OL.newInputDatasetBuilder()
            .namespace(path.toFile().parent).name(path.toFile().name).build()) }
        return inputs
    }

    private static List<OpenLineage.InputDataset> generateInputs(List<Path> inputPaths){
        final inputs = new LinkedList<OpenLineage.InputDataset>()
        inputPaths.each{ path -> inputs.add(OL.newInputDatasetBuilder()
            .namespace(path.toFile().parent).name(path.toFile().name).build()) }
        return inputs
    }

    private static List<OpenLineage.OutputDataset> getOutputsFromTask(TaskRun task){
        final outputs = new LinkedList<OpenLineage.OutputDataset>()
        task.getOutputFilesNames().each{
            outputs.add(OL.newOutputDatasetBuilder()
                .namespace(new URI(task.workDirStr).toString()).name(it).build()) }
        return outputs
    }

    private static List<OpenLineage.OutputDataset> generateOutputs(List<Path> outputsPaths){
        final outputs = new LinkedList<OpenLineage.OutputDataset>()
        outputsPaths.each{ path -> outputs.add(OL.newOutputDatasetBuilder()
                .namespace(path.getParent().toUri().toString())
                .name(path.getName())
                .build())
        }
        return outputs
    }
}