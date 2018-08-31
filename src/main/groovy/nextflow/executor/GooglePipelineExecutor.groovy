package nextflow.executor

import com.google.api.client.googleapis.auth.oauth2.GoogleCredential
import com.google.api.client.googleapis.javanet.GoogleNetHttpTransport
import com.google.api.client.json.jackson2.JacksonFactory
import com.google.api.services.genomics.v2alpha1.Genomics
import com.google.api.services.genomics.v2alpha1.model.*
import com.google.cloud.storage.contrib.nio.CloudStorageConfiguration
import com.google.cloud.storage.contrib.nio.CloudStorageFileSystem
import com.google.cloud.storage.contrib.nio.CloudStorageFileSystemProvider
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

import java.nio.file.Files
import java.nio.file.LinkOption
import java.nio.file.Path
import java.nio.file.Paths
import java.nio.file.attribute.FileAttribute

import static groovy.json.JsonOutput.toJson

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

        log.debug "[GOOGLE PIPELINE] Finished registration for executor $name"
    }

    @Override
    protected TaskMonitor createTaskMonitor() {
        //TODO: Copied directly from the AWS Batch executor. Need to check if these are the values we want for monitoring tasks running on google pipeline.
        TaskPollingMonitor.create(session, name, 1000, Duration.of('10 sec'))
    }

    @Override
    TaskHandler createTaskHandler(TaskRun task) {
        return new GooglePipelineTaskHandler(task, this, pipelineConfig)
    }

    //TODO: Synchronize error messages to be in line with what is used elsewhere in Nextflow
    GooglePipelineConfiguration validateConfiguration() {

        //Make sure that the workdir is a GCE Bucket
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

        new GooglePipelineConfiguration(session.config.navigate("gce.project") as String, session.config.navigate("gce.zone") as String, session.config.navigate("cloud.instanceType") as String)
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
@CompileStatic
class GooglePipelineFileCopyStrategy extends SimpleFileCopyStrategy {

    GooglePipelineTaskHandler handler

    GooglePipelineFileCopyStrategy(TaskBean bean, GooglePipelineTaskHandler handler) {
        super(bean)
        this.handler = handler
    }

    @Override
    String getStageInputFilesScript(Map<String, Path> inputFiles) {
        handler.inputFiles = inputFiles
        "# Google pipeline staging is done in a container that is run before the main container"
    }

    @Override
    String getUnstageOutputFilesScript(List<String> outputFiles, Path targetDir) {
        handler.outputFiles = outputFiles
        handler.outputTargetDir = targetDir
        ": # Google pipeline unstaging is done in a container that is run after the main container"
    }

    //TODO: Pass back to the handler to create the touch file in a the staging action
    @Override
    String touchFile(Path file) {
        ""
        //return super.touchFile(file)
    }
}


@Slf4j
@CompileStatic
class GooglePipelineHelper {

    static String sanitizeName(String name) {
        name.replaceAll(/[^a-zA-Z0-9\-_]+/, '-').take(63)
    }


    static Genomics createGenomicClient() {

        //TODO: Combine shared code with GoogleCLoudDriver into a generic helper
        def credentials = GoogleCredential.applicationDefault

        if (credentials.createScopedRequired()) {
            credentials =
                    credentials.createScoped(["https://www.googleapis.com/auth/cloud-platform"])
        }

        new Genomics.Builder(GoogleNetHttpTransport.newTrustedTransport(), JacksonFactory.defaultInstance, credentials).setApplicationName("NextCode-Experiments/0.1").build()
    }

    //TODO: enum all the flags!!!!
    //TODO: input environment
    static Action createAction(String name, String imageUri, List<String> commands, List<Mount> mounts, List<String> flags = []) {
        new Action().setName(name).setImageUri(imageUri).setCommands(commands).setMounts(mounts).setFlags(flags).setTimeout("3600s")
    }

    static Pipeline createPipeline(List<Action> actions, Resources resources) {
        new Pipeline().setActions(actions).setResources(resources)
    }

    //TODO: constanti-nize the scope and take make it an input.
    static Resources configureResources(String instanceType, String projectId, String zone, String diskName) {

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

/**
 * Implements BASH launcher script for Google Pipeline
 */
//TODO: This code is copied nearly 100% from AWS batch.  Need to grok it and rewrite it
@CompileStatic
class GooglePipelineScriptLauncher extends BashWrapperBuilder {

    GooglePipelineScriptLauncher(TaskBean bean, GooglePipelineTaskHandler handler) {
        super(bean, new GooglePipelineFileCopyStrategy(bean, handler))

        //change the targetDir to be inside of the container
        bean.targetDir = FileHelper.asPath("/work")

        // enable the copying of output file to the GS work dir
        //scratch = true
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
//TODO: Standardize loggin
//TODO: Refaaaactor
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
    final static String fileCopyImage = "google/cloud-sdk:alpine"

    Mount sharedMount
    Pipeline taskPipeline

    @PackageScope
    private Operation operation
    Map<String, Path> inputFiles
    List<String> outputFiles
    Path outputTargetDir
    Metadata metadata


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

        //TODO: These patches should be removed once various kinks in the CloudStorageFilesystemProvider have been solved
        CloudStoragePatcher.patch()

        log.debug "[GOOGLE PIPELINE] Created handler for task '${task.name}'."
    }

    void validateConfiguration() {
        if (!task.container) {
            throw new ProcessUnrecoverableException("No container is specified for process $task.name . Either specify the container to use in the process definition or with 'process.container' value in your config")
        }
    }

    @Override
    boolean checkIfRunning() {

        //log.debug "[GOOGLE PIPELINE] Checking if task '$task.name' is still running"

        operation = executor.genomicsClient.projects().operations().get(operation.getName()).execute()

        //log.debug "[GOOGLE PIPELINE] Task '$task.name' still running = ${!operation.getDone()}"

        return !operation.getDone()
    }

    @Override
    boolean checkIfCompleted() {

        //log.debug "[GOOGLE PIPELINE] Check if task '$task.name' has completed"

        operation = executor.genomicsClient.projects().operations().get(operation.getName()).execute()

        def events = extractRuntimeDataFromOperation()
        events.reverse().each {
            log.trace "[GOOGLE PIPELINE] New event for task $task.name - time: ${it.get("timestamp")} - ${it.get("description")}"
        }

        if (operation.getDone()) {
            log.debug "[GOOGLE PIPELINE] Task '$task.name' complete. Start Time: ${metadata.getStartTime()} - End Time: ${metadata.getEndTime()}"
            //log.debug(operation.toPrettyString())
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
        if(!this.metadata) {
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

        //Create the mount for out work files.
        sharedMount = GooglePipelineHelper.configureMount("workDisk", mountPath)

        def mainAction = GooglePipelineHelper.createAction(taskInstanceName, task.container, ["bash", "-c", "cd " + task.workDir.toString() + " ; ./" + TaskRun.CMD_RUN], [sharedMount])

        def resources = GooglePipelineHelper.configureResources(pipeConfig.vmInstanceType, pipeConfig.project, pipeConfig.zone, "workDisk")

        def stagingCopy = inputFiles.collect {
            "gsutil cp -P ${it.value.isDirectory() ? "-r" : ""} ${it.value.toUriString()} ${task.workDir.toString()}$it.key${it.value.isDirectory() ? "/" : ""} || true"
        }

        //TODO: chmod 777 is bad m'kay
        def stagingScript = """
           mkdir -p $task.workDir ;
           chmod 777 $task.workDir ;
           ${stagingCopy.join(" ; ")} ;
           cd ${task.workDir.toString()} ;
           chmod +x .command*           
        """

        def stagingAction = GooglePipelineHelper.createAction("$taskInstanceName-staging".toString(), fileCopyImage, ["bash", "-c", stagingScript.toString()], [sharedMount], ["ALWAYS_RUN", "IGNORE_EXIT_STATUS"])

        def unstagingCopy = outputFiles.collect {
            "gsutil cp -P -r ${task.workDir}$it ${task.workDir.toUriString()} || true"
        }


        //TODO: Only copy the google directory if we're in deep debug mode
        //TODO: See if we can't skip using blanket copies of everything back to the working directory
        def unstagingScript = """
            ${unstagingCopy.join(" ; ")}
            gsutil cp -P ${task.workDir.toString()}.command* ${task.workDir.toUriString().toString()} || true; 
            gsutil cp -P ${task.workDir.toString()}.exitcode ${task.workDir.toUriString().toString()} || true;
            gsutil cp -P -r /google/ ${task.workDir.toUriString().toString()} || true ;
            gsutil cp -P /google/action/2/stdout  ${task.workDir.toUriString().toString()} || true

        """

        def unstagingAction = GooglePipelineHelper.createAction("$taskInstanceName-unstaging".toString(), fileCopyImage, ["bash", "-c", unstagingScript.toString()], [sharedMount], ["ALWAYS_RUN", "IGNORE_EXIT_STATUS"])


        //TODO: This doesn't work for some reason
        //def listMount = GooglePipelineHelper.createAction("$taskInstanceName-list".toString(), "ubuntu", ["ls", "-Rhalt", "/work"], [sharedMount], ["ALWAYS_RUN", "IGNORE_EXIT_STATUS"])

        taskPipeline = GooglePipelineHelper.createPipeline([stagingAction, mainAction, unstagingAction], resources)

        //log.debug(toJson(taskPipeline))

        operation = executor.genomicsClient.pipelines().run(new RunPipelineRequest().setPipeline(taskPipeline)).execute()

        //log.debug(operation.toPrettyString())

        log.debug "[GOOGLE PIPELINE] Submitted task '$task.name. Assigned Pipeline operation name = '${operation.getName()}'"
    }
}

@Slf4j
class CloudStoragePatcher {

    def static patch() {
        def oldExists = Files.metaClass.getStaticMetaMethod("exists", Path)
        def oldIsDirectory = Files.metaClass.getStaticMetaMethod("isDirectory", Path)
        def oldRelativize = Path.metaClass.getMetaMethod("relativize", Path)

        Path.metaClass.relativize = { Path path ->
            log.debug "Intercepting relativize"
            if (path in CloudStoragePath && !path.startsWith(delegate.toString()[0])) {
                CloudStorageFileSystem noPseudoDirFs = CloudStorageFileSystem.forBucket(path.toUri().getHost(), CloudStorageConfiguration.builder().usePseudoDirectories(false).build())
                def newPath = noPseudoDirFs.getPath("/" + path.toString())
                oldRelativize.invoke(delegate, newPath)
            } else
                oldRelativize.invoke(delegate, path)

        }

        Files.metaClass.static.createDirectories = { Path dir, FileAttribute<?>... attrs ->
            log.debug "CloudStortageFilesystem has to create the following directory: $dir"
        }

        Files.metaClass.static.exists = { Path path, LinkOption... options ->
            def oldret = oldExists.invoke(null, path, options)
            if (!oldret && path in CloudStoragePath && !path.toString().endsWith("/")) {
                CloudStorageFileSystem noPseudoDirFs = CloudStorageFileSystem.forBucket(path.toUri().getHost(), CloudStorageConfiguration.builder().usePseudoDirectories(false).build())
                def pathStr = path.toString().startsWith("/") ? path.toString().drop(1) : path.toString()
                def newPath = noPseudoDirFs.getPath(pathStr + "/")
                def newRet = oldExists.invoke(null, newPath, options)
                newRet
            } else
                oldret
        }

        Files.metaClass.static.isDirectory = { Path path, LinkOption... options ->
            def oldRet = oldIsDirectory.invoke(null, path, options)

            if (!oldRet && path in CloudStoragePath) {
                CloudStorageFileSystem noPseudoDirFs = CloudStorageFileSystem.forBucket(path.toUri().getHost(), CloudStorageConfiguration.builder().usePseudoDirectories(false).build())
                def pathStr = path.toString()
                def newPath = noPseudoDirFs.getPath(!pathStr.endsWith("/") ? pathStr + "/" : pathStr)
                def newRet = oldExists.invoke(null, newPath, options)
                newRet
            } else
                oldRet
        }

        log.debug "Finished patching"
    }
}

