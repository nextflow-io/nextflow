package nextflow.prov.openlineage


import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import io.openlineage.client.OpenLineage
import io.openlineage.client.OpenLineageClient
import io.openlineage.client.transports.HttpTransport
import nextflow.Session
import nextflow.processor.TaskHandler
import nextflow.processor.TaskRun
import nextflow.trace.TraceObserver
import nextflow.trace.TraceRecord

import java.nio.file.Path
import java.time.ZonedDateTime

/**
 * OpenLineage events observer
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Slf4j
@CompileStatic
class OpenLineageObserver implements TraceObserver {
    private OpenLineageClient client
    private OpenLineageEventFactory factory
    private Session session

    OpenLineageObserver(String server){
        this.client = OpenLineageClient.builder()
                .transport(HttpTransport.builder().uri(server).build())
                .build();
    }

    @Override
    void onFlowCreate(Session session) {
        this.session = session
        this.factory = new OpenLineageEventFactory()
        this.factory.setWorkflowJobFromSession(session)
        this.factory.setWorkflowRun(session.uniqueId)
    }

    @Override
    void onFlowComplete() {
        if( factory ) {
            final event = factory.createWorkflowCompleteEvent()
            client.emit(event)
            log.debug("[Openlineage] Workflow complete event sent.")
        }
    }

    @Override
    void onFlowBegin(){
        final event = factory.createWorkflowStartEvent()
        client.emit(event)
        log.debug("[Openlineage] Workflow start event sent.")
    }

    @Override
    void onProcessStart(TaskHandler handler, TraceRecord trace){
        final event = factory.createTaskStartEvent(handler.task, trace)
        client.emit(event)
        log.debug("[Openlineage] Task ${handler.task.id} start event sent.")
    }

    @Override
    void onProcessComplete(TaskHandler handler, TraceRecord trace){
        final event = factory.createTaskCompleteEvent(handler.task, trace)
        client.emit(event)
        log.debug("[Openlineage] Task ${handler.task.id} complete event sent.")
    }

    @Override
    void onProcessCached(TaskHandler handler, TraceRecord trace){
        final event = factory.createTaskCachedEvent(handler.task, trace)
        client.emit(event)
        log.debug("[Openlineage] Task ${handler.task.id} cached event sent.")
    }

    @Override
    void onFlowError(TaskHandler var1, TraceRecord var2){
        if (factory) {
            final event = factory.createWorkflowFailedEvent()
            client.emit(event)
            log.debug("[Openlineage] Workflow error event sent.")
        }
    }

    @Override
    void onFilePublish(Path destination, Path source){
        final event = factory.createPublishEvent(source, destination)
        client.emit(event)
        log.debug("[Openlineage] File publish ")
    }
}