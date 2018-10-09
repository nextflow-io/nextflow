package nextflow.executor

import com.google.api.client.googleapis.auth.oauth2.GoogleCredential
import com.google.api.client.googleapis.javanet.GoogleNetHttpTransport
import com.google.api.client.json.jackson2.JacksonFactory
import com.google.api.services.genomics.v2alpha1.Genomics
import com.google.api.services.genomics.v2alpha1.model.*
import com.google.cloud.storage.contrib.nio.CloudStoragePath
import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.exception.AbortOperationException
import nextflow.exception.ProcessUnrecoverableException
import nextflow.file.FileHelper
import nextflow.processor.*
import nextflow.script.ScriptType
import nextflow.util.Duration

import java.nio.file.Path

@Slf4j
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
        log.debug "[GOOGLE PIPELINE] Pipeline config: $pipelineConfig"

        genomicsClient = GooglePipelineHelper.createGenomicClient()

        log.debug "[GOOGLE PIPELINE] Finished registration for executor $name"
    }

    @Override
    protected TaskMonitor createTaskMonitor() {
        TaskPollingMonitor.create(session, name, 1000, Duration.of('10 sec'))
    }

    @Override
    TaskHandler createTaskHandler(TaskRun task) {
        return new GooglePipelineTaskHandler(task, this, pipelineConfig)
    }

    GooglePipelineConfiguration validateConfiguration() {

        //Make sure that the workdir is a GS Bucket
        if (!(session.workDir instanceof CloudStoragePath)) {
            session.abort()
            throw new AbortOperationException("When using `$name` executor a GCE bucket must be provided as a working directory -- Add the option `-w gs://<your-bucket/path>` to your run command line or specify a workDir in your config file.")
        }

        //Check for the existence of all required configuration for our executor
        def requiredConfig = ["gce.project", "gce.zone"]

        requiredConfig.each {
            if (!session.config.navigate(it)) {
                session.abort()
                throw new AbortOperationException("Required config value '$it' for executor $name is not defined. Please add it to your process or nextflow configuration file.")
            }
        }

        new GooglePipelineConfiguration(
                session.config.navigate("gce.project") as String,
                session.config.navigate("gce.zone") as String,
                session.config.navigate("cloud.instanceType") as String,
                session.config.navigate("cloud.preemptible") as boolean
        )
    }
}


@CompileStatic
class GooglePipelineConfiguration {
    String project
    String zone
    String vmInstanceType
    boolean preemptible

    GooglePipelineConfiguration(String project, String zone, String vmInstanceType, boolean preemptible = false) {
        this.project = project
        this.zone = zone
        this.vmInstanceType = vmInstanceType
        this.preemptible = preemptible
    }


    @Override
    String toString() {
        "GooglePipelineConfiguration{" +
                "project='" + project + '\'' +
                ", zone='" + zone + '\'' +
                ", vmInstanceType='" + vmInstanceType + '\'' +
                ", preemptible=" + preemptible +
                '}'
    }
}

@Slf4j
@CompileStatic
class GooglePipelineFileCopyStrategy extends SimpleFileCopyStrategy {

    GooglePipelineTaskHandler handler
    TaskBean task

    GooglePipelineFileCopyStrategy(TaskBean bean, GooglePipelineTaskHandler handler) {
        super(bean)
        this.handler = handler
        this.task = bean
    }

    @Override
    String getStageInputFilesScript(Map<String, Path> inputFiles) {

        def stagingCommands = inputFiles.collect {
            "gsutil -m  -q cp -P -c -r ${it.value.toUriString()} ${task.workDir}/${it.value.toString().endsWith(it.key) ? "" : it.key} || true".toString()
        }

        log.debug "[GOOGLE PIPELINE] Constructed the following file copy staging commands: $stagingCommands"

        handler.stagingCommands.addAll(stagingCommands)

        //Insert this comment into the task run script to note that the staging is done differently
        "# Google pipeline staging is done in a pipeline action step that is run prior to the main pipeline action"
    }

    @Override
    String getUnstageOutputFilesScript(List<String> outputFiles, Path targetDir) {

        def unstagingCommands = outputFiles.collect {
            "gsutil -m -q cp -P -r -c ${task.workDir}/$it ${task.workDir.toUriString()} || true".toString()
        }

        log.debug "[GOOGLE PIPELINE] Constructed the following file copy staging commands: $unstagingCommands"

        handler.unstagingCommands.addAll(unstagingCommands)

        //Insert this comment into the task run script to note that the unstaging is done differently
        ": # Google pipeline unstaging is done in a pipeline action step that is run after the main pipeline action"
    }

    //Although it seems like this is used as a general "touch" mechanism it's actually just used to create a .begin in the working directory
    @Override
    String touchFile(Path file) {
        if (file in CloudStoragePath) {
            handler.stagingCommands << "echo start | gsutil -q cp  -c - ${file.toUriString()} || true".toString()
            "# Google pipeline touchFile is done in a pipeline action step that is run prior to the main pipeline action"
        } else
            return super.touchFile(file)
    }
}

@Slf4j
@CompileStatic
class GooglePipelineHelper {

    static final String SCOPE_CLOUD_PLATFORM = "https://www.googleapis.com/auth/cloud-platform"
    static final List<String> ENV_VAR_TO_INCLUDE = ["NXF_DEBUG"]

    enum ActionFlags {
        FLAG_UNSPECIFIED,
        IGNORE_EXIT_STATUS,
        RUN_IN_BACKGROUND,
        ALWAYS_RUN,
        ENABLE_FUSE,
        PUBLISH_EXPOSED_PORTS,
        DISABLE_IMAGE_PREFETCH,
        DISABLE_STANDARD_ERROR_CAPTURE
    }

    static String sanitizeName(String name) {
        name.replaceAll(/[^a-zA-Z0-9\-_]+/, '-').take(63)
    }

    static Genomics createGenomicClient() {

        def credentials = GoogleCredential.applicationDefault

        if (credentials.createScopedRequired()) {
            credentials = credentials.createScoped([SCOPE_CLOUD_PLATFORM])
        }

        new Genomics.Builder(GoogleNetHttpTransport.newTrustedTransport(), JacksonFactory.defaultInstance, credentials)
                .setApplicationName("Nextflow GooglePipelineExecutor")
                .build()
    }

    static Map<String,String> getEnvironment() {
        def ret = [:]
        def env = System.getenv()

        env.each { kv ->
            if(ENV_VAR_TO_INCLUDE.contains(kv.key))
                ret << kv
        }
        ret
    }

    static Action createAction(String name, String imageUri, List<String> commands, List<Mount> mounts, List<ActionFlags> flags = [], String entrypoint = null) {
        new Action()
                .setName(name)
                .setImageUri(imageUri)
                .setCommands(commands)
                .setMounts(mounts)
                .setFlags(flags.collect{flag -> flag.toString()})
                .setEntrypoint(entrypoint)
                .setEnvironment(getEnvironment())
    }

    static Pipeline createPipeline(List<Action> actions, Resources resources) {
        new Pipeline().setActions(actions).setResources(resources)
    }

    //TODO: Do we want to configure this via nextflow config?
    static Resources configureResources(String instanceType, String projectId, String zone, String diskName, List<String> scopes = null,boolean preEmptible = false) {

        def disk = new Disk()
        disk.setName(diskName)

        def serviceAccount = new ServiceAccount()
        if (scopes)
            serviceAccount.setScopes(scopes)

        def vm = new VirtualMachine()
                .setMachineType(instanceType)
                .setDisks([disk])
                .setServiceAccount(serviceAccount)
                .setPreemptible(preEmptible)


        new Resources()
                .setProjectId(projectId)
                .setZones([zone])
                .setVirtualMachine(vm)
    }

    static Mount configureMount(String diskName, String mountPath, boolean readOnly = false) {
        new Mount().setDisk(diskName).setPath(mountPath).setReadOnly(readOnly)
    }

}

/**
 * Implements BASH launcher script for Google Pipeline
 */
//TODO: This code is copied nearly 100% from AWS batch.  Need to grok it and rewrite it
//TODO: Scratch stuff and changing the targetDir are both unsuitable to trigger getUnstageOutputFilesScript
@CompileStatic
class GooglePipelineScriptLauncher extends BashWrapperBuilder {

    GooglePipelineScriptLauncher(TaskBean bean, GooglePipelineTaskHandler handler) {
        super(bean, new GooglePipelineFileCopyStrategy(bean, handler))

        // enable the copying of output file to the GS work dir
        //scratch = false

        //TODO: here lies a hackdragon, 'YARR!!!!
        bean.targetDir = FileHelper.asPath("/work")

        // include task script as an input to force its staging in the container work directory
        bean.inputFiles[TaskRun.CMD_SCRIPT] = bean.workDir.resolve(TaskRun.CMD_SCRIPT)
        // include the wrapper script as in input to force its staging in the container work directory
        bean.inputFiles[TaskRun.CMD_RUN] = bean.workDir.resolve(TaskRun.CMD_RUN)
        // add the wrapper file when stats are enabled
        if (bean.statsEnabled) {
            bean.inputFiles[TaskRun.CMD_STUB] = bean.workDir.resolve(TaskRun.CMD_STUB)
        }
        // include task stdin file
        if (bean.input != null) {
            bean.inputFiles[TaskRun.CMD_INFILE] = bean.workDir.resolve(TaskRun.CMD_INFILE)
        }
    }

}


@Slf4j
//@CompileStatic
class GooglePipelineTaskHandler extends TaskHandler {

    final GooglePipelineExecutor executor
    final TaskBean taskBean
    final GooglePipelineConfiguration pipeConfig

    final String taskName
    final String taskInstanceName

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
    final static String diskName = "nf-pipeline-work"
    final static String fileCopyImage = "google/cloud-sdk:alpine"

    Mount sharedMount
    Pipeline taskPipeline

    private Operation operation
    private Metadata metadata

    @PackageScope
    List<String> stagingCommands = []
    @PackageScope
    List<String> unstagingCommands = []

    GooglePipelineTaskHandler(TaskRun task, Executor executor, GooglePipelineConfiguration pipeConfig) {
        super(task)
        this.executor = executor as GooglePipelineExecutor
        this.taskBean = new TaskBean(task)
        this.pipeConfig = pipeConfig

        this.logFile = task.workDir.resolve(TaskRun.CMD_LOG)
        this.scriptFile = task.workDir.resolve(TaskRun.CMD_SCRIPT)
        this.inputFile = task.workDir.resolve(TaskRun.CMD_INFILE)
        this.outputFile = task.workDir.resolve(TaskRun.CMD_OUTFILE)
        this.errorFile = task.workDir.resolve(TaskRun.CMD_ERRFILE)
        this.exitFile = task.workDir.resolve(TaskRun.CMD_EXIT)
        this.wrapperFile = task.workDir.resolve(TaskRun.CMD_RUN)
        this.stubFile = task.workDir.resolve(TaskRun.CMD_STUB)
        this.traceFile = task.workDir.resolve(TaskRun.CMD_TRACE)

        this.taskName = GooglePipelineHelper.sanitizeName("nf-task-${executor.session.uniqueId}-${task.name}")
        this.taskInstanceName = GooglePipelineHelper.sanitizeName("$taskName-$task.id")

        validateConfiguration()

        log.debug "[GOOGLE PIPELINE] Created handler for task '${task.name}'."
    }

    void validateConfiguration() {
        if (!task.container) {
            throw new ProcessUnrecoverableException("No container is specified for process $task.name . Either specify the container to use in the process definition or with 'process.container' value in your config")
        }
    }

    @Override
    boolean checkIfRunning() {
        operation = executor.genomicsClient.projects().operations().get(operation.getName()).execute()
        return !operation.getDone()
    }

    @Override
    //TODO: Catch pipeline errors and report them back
    boolean checkIfCompleted() {
        operation = executor.genomicsClient.projects().operations().get(operation.getName()).execute()

        def events = extractRuntimeDataFromOperation()
        events.reverse().each {
            log.trace "[GOOGLE PIPELINE] New event for task '$task.name' - time: ${it.get("timestamp")} - ${it.get("description")}"
        }

        if (operation.getDone()) {
            log.debug "[GOOGLE PIPELINE] Task '$task.name' complete. Start Time: ${metadata.getStartTime()} - End Time: ${metadata.getEndTime()}"

            // finalize the task
            task.exitStatus = readExitFile()
            task.stdout = outputFile
            task.stderr = errorFile
            status = TaskStatus.COMPLETED
            return true
        } else
            return false
    }

    private List<Event> extractRuntimeDataFromOperation() {
        def metadata = operation.getMetadata() as Metadata
        if (!this.metadata) {
            this.metadata = metadata
            return metadata.getEvents()
        } else {
            //Get the new events
            def delta = metadata.getEvents().size() - this.metadata.getEvents().size()
            this.metadata = metadata
            return delta > 0 ? metadata.getEvents().take(delta) : [] as List<Event>
        }
    }

    private int readExitFile() {
        try {
            exitFile.text as Integer
        }
        catch (Exception e) {
            log.debug "[GOOGLE PIPELINE] Cannot read exitstatus for task: `$task.name`", e
            return Integer.MAX_VALUE
        }
    }

    @Override
    void kill() {
        log.debug "[GOOGLE PIPELINE] Killing pipeline '$operation.name'"
        executor.genomicsClient.projects().operations().cancel(operation.getName(), new CancelOperationRequest())
    }

    @Override
    void submit() {

        final launcher = new GooglePipelineScriptLauncher(this.taskBean, this)
        launcher.build()

        //TODO: chmod 777 is bad m'kay
        //TODO: eliminate cd command as well as wildcard +x
        def stagingScript = """
           mkdir -p $task.workDir ;
           chmod 777 $task.workDir ;
           ${stagingCommands.join(" ; ")} ;
           cd $task.workDir ;
           chmod 777 ${TaskRun.CMD_SCRIPT} ${TaskRun.CMD_RUN}            
        """.stripIndent().leftTrim()

        def mainScript = "cd ${task.workDir} ; echo \$(./${TaskRun.CMD_RUN}) | bash 2>&1 | tee ${TaskRun.CMD_LOG}".toString()

        /*
         * -m = run in parallel
         * -q = quiet mode
         * cp = copy
         * -P = preserve POSIX attributes
         * -c = continues on errors
         * -r = recursive copy
         */
        def gsCopyPrefix = "gsutil -m -q cp -P -c"

        //Copy the logs provided by Google Pipelines for the pipline to our work dir.
        if(System.getenv().get("NXF_DEBUG")) {
            unstagingCommands << "$gsCopyPrefix -r /google/ ${task.workDir.toUriString()} || true".toString()
        }

        //add the task output files to unstaging command list
        [TaskRun.CMD_ERRFILE,
         TaskRun.CMD_OUTFILE,
         TaskRun.CMD_EXIT,
         TaskRun.CMD_LOG
        ].each {
            unstagingCommands << "$gsCopyPrefix ${task.workDir}/$it ${task.workDir.toUriString()} || true".toString()
        }

        //Copy nextflow task progress files as well as the files we need to unstage
        def unstagingScript = """                                                
            ${unstagingCommands.join(" ; ")}                        
        """.stripIndent().leftTrim()

        log.debug "Staging script for task $task.name -> $stagingScript"
        log.debug "Main script for task $task.name -> $mainScript"
        log.debug "Unstaging script for task $task.name -> $unstagingScript"

        //Create the mount for out work files.
        sharedMount = GooglePipelineHelper.configureMount(diskName, mountPath)

        //need the cloud-platform scope so that we can execute gsutil cp commands
        def resources = GooglePipelineHelper.configureResources(pipeConfig.vmInstanceType, pipeConfig.project, pipeConfig.zone, diskName, [GooglePipelineHelper.SCOPE_CLOUD_PLATFORM],pipeConfig.preemptible)

        def stagingAction = GooglePipelineHelper.createAction("$taskInstanceName-staging", fileCopyImage, ["bash", "-c", stagingScript], [sharedMount], [GooglePipelineHelper.ActionFlags.ALWAYS_RUN, GooglePipelineHelper.ActionFlags.IGNORE_EXIT_STATUS])
        //TODO: Do we really want to override the entrypoint?
        def mainAction = GooglePipelineHelper.createAction(taskInstanceName, task.container, ['-o', 'pipefail', '-c', mainScript], [sharedMount], [], "bash")

        def unstagingAction = GooglePipelineHelper.createAction("$taskInstanceName-unstaging", fileCopyImage, ["bash", "-c", unstagingScript], [sharedMount], [GooglePipelineHelper.ActionFlags.ALWAYS_RUN, GooglePipelineHelper.ActionFlags.IGNORE_EXIT_STATUS])

        taskPipeline = GooglePipelineHelper.createPipeline([stagingAction, mainAction, unstagingAction], resources)

        operation = executor.genomicsClient.pipelines().run(new RunPipelineRequest().setPipeline(taskPipeline)).execute()

        log.trace "[GOOGLE PIPELINE] Submitted task '$task.name. Assigned Pipeline operation name = '${operation.getName()}'"
    }
}