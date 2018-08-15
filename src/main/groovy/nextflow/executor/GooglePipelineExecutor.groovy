package nextflow.executor

import com.google.api.client.googleapis.auth.oauth2.GoogleCredential
import com.google.api.client.googleapis.javanet.GoogleNetHttpTransport
import com.google.api.client.json.jackson2.JacksonFactory
import com.google.api.services.genomics.v2alpha1.Genomics
import com.google.api.services.genomics.v2alpha1.model.Action
import com.google.api.services.genomics.v2alpha1.model.Disk
import com.google.api.services.genomics.v2alpha1.model.Mount
import com.google.api.services.genomics.v2alpha1.model.Operation
import com.google.api.services.genomics.v2alpha1.model.Pipeline
import com.google.api.services.genomics.v2alpha1.model.Resources
import com.google.api.services.genomics.v2alpha1.model.RunPipelineRequest
import com.google.api.services.genomics.v2alpha1.model.ServiceAccount
import com.google.api.services.genomics.v2alpha1.model.VirtualMachine
import com.google.cloud.storage.contrib.nio.CloudStoragePath
import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.exception.AbortOperationException
import nextflow.exception.ProcessUnrecoverableException
import nextflow.processor.TaskBean
import nextflow.processor.TaskHandler
import nextflow.processor.TaskMonitor
import nextflow.processor.TaskPollingMonitor
import nextflow.processor.TaskRun
import nextflow.script.ScriptType
import nextflow.util.Duration

import java.nio.file.Path

import static groovy.json.JsonOutput.*

@Slf4j
@CompileStatic
@SupportedScriptTypes(ScriptType.SCRIPTLET)
class GooglePipelineExecutor extends Executor {

    @PackageScope
    static Genomics genomicsClient

    static GooglePipelineConfiguration pipelineConfig

    @Override
    final boolean isContainerNative() {
        return true
    }

    @Override
    void register() {
        super.register()

        pipelineConfig = validateConfiguration()

        genomicsClient = GooglePipelineHelper.createGenomicClient()

        log.debug "Finished registration for executor $name"
    }

    @Override
    protected TaskMonitor createTaskMonitor() {
        //TODO: Copied directly from the AWS Batch executor. Need to check if these are the values we want for monitoring tasks running on google pipeline.
        TaskPollingMonitor.create(session, name, 1000, Duration.of('20 sec'))
    }

    @Override
    TaskHandler createTaskHandler(TaskRun task) {
        return new GooglePipelineTaskHandler(task, this,pipelineConfig)
    }

    //TODO: Synchronize error messages to be in line with what is used elsewhere in Nextflow
    GooglePipelineConfiguration validateConfiguration() {

        //Make sure that the workdir is a GCE Bucket
        if (!(session.workDir instanceof CloudStoragePath)) {
            session.abort()
            throw new AbortOperationException("When using `$name` executor a GCE bucket must be provided as a working directory -- Add the option `-w gs://<your-bucket/path>` to your run command line or specify a workDir in your config file.")
        }


        //Check for the existence of all required configuration for our executor
        def requiredConfig = ["gce.project","gce.zone"]

        requiredConfig.each {
            if(!session.config.navigate(it)) {
                session.abort()
                throw new AbortOperationException("Required config value '$it' for executor $name is not defined. Please add it to your process or nextflow configuration file.")
            }
        }

        new GooglePipelineConfiguration(session.config.navigate("gce.project") as String,session.config.navigate("gce.zone") as String,session.config.navigate("cloud.instanceType") as String)
    }
}


class GooglePipelineConfiguration {
    String project
    String zone
    String vmInstanceType

    GooglePipelineConfiguration(String project, String zone, String vmInstanceType) {
        this.project = project
        this.zone = zone
        this.vmInstanceType = vmInstanceType
    }
}


@Slf4j
class GooglePipelineHelper {


    static Genomics createGenomicClient() {

        //TODO: Do we want a different way to get the credentials?
        //TODO: Combine shared code with GoogleCLoudDriver into a generic helper
        def credentials = GoogleCredential.applicationDefault

        if (credentials.createScopedRequired()) {
            credentials =
                    credentials.createScoped(["https://www.googleapis.com/auth/cloud-platform"])
        }

        new Genomics.Builder(GoogleNetHttpTransport.newTrustedTransport(), JacksonFactory.defaultInstance, credentials).setApplicationName("NextCode-Experiments/0.1").build()
    }


    //TODO: enum all the flags!!!!
    static Action createAction(String name, String imageUri, List<String> commands, List<Mount> mounts, List<String> flags = []) {
        new Action().setName(name).setImageUri(imageUri).setCommands(commands).setMounts(mounts).setFlags(flags)
    }

    static Pipeline createPipeline(List<Action> actions,Resources resources) {
        new Pipeline().setActions(actions).setResources(resources)
    }


    //TODO: constanti-nize the scope and take make it an input.
    static Resources configureResources(String instanceType,String projectId,String zone,String diskName) {

        def disk = new Disk()
            disk.setName(diskName)


        def serviceAccount = new ServiceAccount()
            .setScopes(["https://www.googleapis.com/auth/cloud-platform"])


        def vm = new VirtualMachine()
            .setMachineType(instanceType)
            .setDisks([disk])
            .setServiceAccount(serviceAccount)



        new Resources()
            .setProjectId(projectId)
            .setZones([zone])
            .setVirtualMachine(vm)
    }

    static Mount configureMount(String diskName, String mountPath, boolean readOnly = false) {
        new Mount().setDisk(diskName).setPath(mountPath).setReadOnly(readOnly)
    }

}

@Slf4j
//@CompileStatic
//TODO: Enable static again and reeefactor
class GooglePipelineTaskHandler extends TaskHandler {

    final GooglePipelineExecutor executor
    final TaskBean taskBean
    final GooglePipelineConfiguration pipeConfig


    private final Path exitFile
    private final Path wrapperFile
    private final Path outputFile
    private final Path errorFile
    private final Path logFile
    private final Path scriptFile
    private final Path inputFile
    private final Path stubFile
    private final Path traceFile

    //Constants
    final static String mountPath = "/work"
    final static String fileCopyImage = "google/cloud-sdk"

    Mount sharedMount
    Pipeline taskPipeline

    @PackageScope
    private Operation operation


    GooglePipelineTaskHandler(TaskRun task, Executor executor,GooglePipelineConfiguration pipeConfig) {
        super(task)
        this.executor = executor as GooglePipelineExecutor
        this.taskBean = new TaskBean(task)
        this.pipeConfig = pipeConfig

        this.logFile = task.workDir.resolve(TaskRun.CMD_LOG)
        this.scriptFile = task.workDir.resolve(TaskRun.CMD_SCRIPT)
        this.inputFile =  task.workDir.resolve(TaskRun.CMD_INFILE)
        this.outputFile = task.workDir.resolve(TaskRun.CMD_OUTFILE)
        this.errorFile = task.workDir.resolve(TaskRun.CMD_ERRFILE)
        this.exitFile = task.workDir.resolve(TaskRun.CMD_EXIT)
        this.wrapperFile = task.workDir.resolve(TaskRun.CMD_RUN)
        this.stubFile = task.workDir.resolve(TaskRun.CMD_STUB)
        this.traceFile = task.workDir.resolve(TaskRun.CMD_TRACE)

        validateConfiguration()


        sharedMount = GooglePipelineHelper.configureMount(generateTaskName(),mountPath)

        def mainAction = GooglePipelineHelper.createAction(generateTaskInstanceName(),task.container,["bash","-c","echo liverpool ; exit 4"],[sharedMount])

        def resources = GooglePipelineHelper.configureResources(pipeConfig.vmInstanceType,pipeConfig.project,pipeConfig.zone,generateTaskName())


        def postCopyFiles = [ "/google/logs/action/1/stdout": logFile.toUriString(),
                              "/google/logs/action/1/stderr": errorFile.toUriString()
        ]


        def postCopyAction = genererateFileCopyAction("post-copy",postCopyFiles)



        taskPipeline = GooglePipelineHelper.createPipeline([mainAction,postCopyAction],resources)

        log.debug "Created task ${task.name}."
    }



    Action genererateFileCopyAction(String nameSuffix,  Map<String,String> copies) {

        List<String> commands = copies.collect {
            "gsutil cp $it.key $it.value".toString()
        }

        String comb = commands.join(" ; ")


        GooglePipelineHelper.createAction("${generateTaskInstanceName()}-$nameSuffix",fileCopyImage,["bash","-c",comb].toList(),[sharedMount],["ALWAYS_RUN","IGNORE_EXIT_STATUS"])

    }



    String generateTaskName() {
        "nf-task-${executor.session.uniqueId}-${task.name}"
    }


    String generateTaskInstanceName() {
        "${generateTaskName()}-${task.id}"
    }

    void validateConfiguration() {
        if(!task.container) {
            throw new ProcessUnrecoverableException("No container is specified for process $task.name . Either specify the container to use in the process definition or with 'process.container' value in your config")
        }
    }
www
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


        if (operation.getDone()) {
            log.debug(operation.toPrettyString())
            if (operation.getError()) {
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
