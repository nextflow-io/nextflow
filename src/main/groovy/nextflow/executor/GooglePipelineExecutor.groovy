package nextflow.executor

import com.google.api.client.googleapis.auth.oauth2.GoogleCredential
import com.google.api.client.googleapis.javanet.GoogleNetHttpTransport
import com.google.api.client.json.jackson2.JacksonFactory
import com.google.api.services.genomics.v2alpha1.Genomics
import com.google.api.services.genomics.v2alpha1.model.Action
import com.google.api.services.genomics.v2alpha1.model.Operation
import com.google.api.services.genomics.v2alpha1.model.Pipeline
import com.google.api.services.genomics.v2alpha1.model.Resources
import com.google.api.services.genomics.v2alpha1.model.RunPipelineRequest
import com.google.api.services.genomics.v2alpha1.model.VirtualMachine
import com.google.cloud.storage.contrib.nio.CloudStoragePath
import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.Global
import nextflow.exception.AbortOperationException
import nextflow.processor.TaskBean
import nextflow.processor.TaskHandler
import nextflow.processor.TaskMonitor
import nextflow.processor.TaskPollingMonitor
import nextflow.processor.TaskRun
import nextflow.util.Duration

import static groovy.json.JsonOutput.*

@Slf4j
@CompileStatic
class GooglePipelineExecutor extends Executor {

    @PackageScope
    static Genomics genomicsClient

    @Override
    void register() {
        super.register()

        //Make sure that the workdir is a GCE Bucket
        if( !(session.workDir instanceof CloudStoragePath) ) {
            session.abort()
            throw new AbortOperationException("When using `$name` executor a GCE bucket must be provided as a working directory -- Add the option `-w gs://<your-bucket/path>` to your run command line")
        }


        def httpTransport = GoogleNetHttpTransport.newTrustedTransport()


        def jsonFactory = JacksonFactory.defaultInstance

        //TODO: Do we want a different way to get the credentials?
        //TODO: Combine shared code with GoogleCLoudDriver into a generic helper
        def credentials = GoogleCredential.applicationDefault

        if (credentials.createScopedRequired()) {
            credentials =
                    credentials.createScoped(["https://www.googleapis.com/auth/cloud-platform"])
        }

        log.debug "Creating google genomics client with following parameters: transport $httpTransport - jsonfactory: $jsonFactory - cred: $credentials"


        def builder = new Genomics.Builder(httpTransport, jsonFactory,credentials).setApplicationName("NextCode-Experiments/0.1")


        genomicsClient = builder.build()

        log.debug "Finished google pipeline registration"
    }

    @Override
    protected TaskMonitor createTaskMonitor() {
        //TODO: Copied directly from the AWS Batch executor. Need to check if these are the values we want for monitoring tasks running on google pipeline.
        TaskPollingMonitor.create(session, name, 1000, Duration.of('20 sec'))
    }

    @Override
    TaskHandler createTaskHandler(TaskRun task) {
        return new GooglePipelineTaskHandler(task,this)
    }
}

@Slf4j
//@CompileStatic
//TODO: Enable static again and reeefactor
class GooglePipelineTaskHandler extends TaskHandler {

    final GooglePipelineExecutor executor
    final TaskBean taskBean

    Pipeline taskPipeline

    private Operation operation


    GooglePipelineTaskHandler(TaskRun task, Executor executor) {
        super(task)
        this.executor = executor as GooglePipelineExecutor
        this.taskBean = new TaskBean(task)

        taskPipeline = new Pipeline()

        def action = new Action()
        action.setName("nf-task-${task.name}")
        action.setImageUri(task.container)
        //Make the container fail for tests
        //action.setCommands(["babayaga!"])
        action.setCommands(["/bin/sleep","1m"])

        taskPipeline.setActions([action])
        taskPipeline.setResources(new Resources().setProjectId(Global.config?.gce?.project).setZones([Global.config?.gce?.zone]).setVirtualMachine(new VirtualMachine().setMachineType(Global.config?.cloud?.instanceType)))

        log.debug "Created task ${taskBean.name}."
    }

    @Override
    boolean checkIfRunning() {

        log.debug "Check if still Running"

        operation = executor.genomicsClient.projects().operations().get(operation.getName()).execute()

        log.debug(toJson(operation.getMetadata().events))

        return !operation.getDone()
    }

    @Override
    boolean checkIfCompleted() {

        log.debug "Check if Completed"

        operation = executor.genomicsClient.projects().operations().get(operation.getName()).execute()


        if(operation.getDone()) {
            log.debug(operation.toPrettyString())
            if(operation.getError()) {
                task.exitStatus = operation.getError().getCode()

            } else {
                task.exitStatus = 0
            }
            return true
        } else {
            log.debug(toJson(operation.getMetadata().events))
            return false
        }

    }

    @Override
    void kill() {

    }

    @Override
    void submit() {

        final launcher = new BashWrapperBuilder(task)
        def p = launcher.build()

        operation = executor.genomicsClient.pipelines().run(new RunPipelineRequest().setPipeline(taskPipeline)).execute()

        log.debug(operation.toPrettyString())

        log.debug "Submit done"
    }
}
