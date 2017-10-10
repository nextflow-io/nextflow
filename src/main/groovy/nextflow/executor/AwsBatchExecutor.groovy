package nextflow.executor
import java.nio.file.Path
import java.nio.file.Paths

import com.amazonaws.services.batch.AWSBatchClient
import com.amazonaws.services.batch.model.CancelJobRequest
import com.amazonaws.services.batch.model.ContainerOverrides
import com.amazonaws.services.batch.model.ContainerProperties
import com.amazonaws.services.batch.model.DescribeJobDefinitionsRequest
import com.amazonaws.services.batch.model.DescribeJobsRequest
import com.amazonaws.services.batch.model.Host
import com.amazonaws.services.batch.model.JobDefinition
import com.amazonaws.services.batch.model.JobDefinitionType
import com.amazonaws.services.batch.model.KeyValuePair
import com.amazonaws.services.batch.model.MountPoint
import com.amazonaws.services.batch.model.RegisterJobDefinitionRequest
import com.amazonaws.services.batch.model.SubmitJobRequest
import com.amazonaws.services.batch.model.Volume
import com.upplication.s3fs.S3Path
import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.Memoized
import groovy.transform.PackageScope
import groovy.transform.ToString
import groovy.util.logging.Slf4j
import nextflow.Nextflow
import nextflow.cloud.aws.AmazonCloudDriver
import nextflow.exception.AbortOperationException
import nextflow.exception.ProcessNotRecoverableException
import nextflow.extension.FilesEx
import nextflow.processor.TaskBean
import nextflow.processor.TaskHandler
import nextflow.processor.TaskMonitor
import nextflow.processor.TaskPollingMonitor
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import nextflow.processor.TaskStatus
import nextflow.util.CacheHelper
import nextflow.util.Duration
import nextflow.util.Escape
/**
 * AWS Batch executor
 * https://aws.amazon.com/batch/
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class AwsBatchExecutor extends Executor {

    /**
     * AWS batch client instance
     */
    @PackageScope
    private static AWSBatchClient client

    private static Path remoteBinDir = null

    final boolean isContainerNative() {
        return true
    }

    @Override
    void register() {
        super.register()

        /*
         * make sure the work dir is a S3 bucket
         */
        if( !(session.workDir instanceof S3Path) ) {
            session.abort()
            throw new AbortOperationException("When using `$name` executor a S3 bucket must be provided as working directory -- Add the option `-w s3://<your-bucket/path>` to your run command line")
        }

        def path = session.config.navigate('env.PATH')
        if( path ) {
            log.warn "Environment PATH defined in config file is ignored by AWS Batch executor"
        }

        /*
         * upload local binaries
         */
        def disableBinDir = session.getExecConfigProp(name, 'disableRemoteBinDir', false)
        if( session.binDir && !session.binDir.empty() && !disableBinDir ) {
            def s3 = Nextflow.tempDir()
            log.info "Uploading local `bin` scripts folder to ${s3.toUriString()}/bin"
            remoteBinDir = FilesEx.copyTo(session.binDir, s3)
        }

        /*
         * retrieve config and credentials and create AWS client
         */
        client = new AmazonCloudDriver(session.config).getBatchClient()

    }

    @PackageScope
    Path getRemoteBinDir() {
        remoteBinDir
    }

    @PackageScope
    AWSBatchClient getClient() {
        client
    }

    @Override
    protected TaskMonitor createTaskMonitor() {
        TaskPollingMonitor.create(session, name, 50, Duration.of('5 sec'))
    }

    @Override
    TaskHandler createTaskHandler(TaskRun task) {
        assert task
        assert task.workDir
        log.trace "[AWS BATCH] Launching process > ${task.name} -- work folder: ${task.workDirStr}"
        new AwsBatchTaskHandler(task, this)
    }
}

/**
 * Implements a task handler for AWS Batch jobs
 */
@Slf4j
@CompileStatic
class AwsBatchTaskHandler extends TaskHandler {

    private final Path exitFile

    private final Path wrapperFile

    private final Path outputFile

    private final Path errorFile

    private final Path logFile

    private final Path scriptFile

    private final Path inputFile

    private final Path stubFile

    private final Path traceFile

    private TaskBean bean

    private AwsBatchExecutor executor

    private AWSBatchClient client

    private volatile String jobId

    private Map<String,String> environment

    final static private Map<String,String> jobDefinitions = [:]


    /** only for testing purpose -- do not use */
    protected AwsBatchTaskHandler() {}

    /**
     * Create a new Batch task handler
     *
     * @param task The {@link TaskRun} descriptor of the task to run
     * @param executor The {@link AwsBatchExecutor} instance
     */
    AwsBatchTaskHandler(TaskRun task, AwsBatchExecutor executor) {
        super(task)
        this.bean = new TaskBean(task)
        this.executor = executor
        this.client = executor.client
        this.environment = System.getenv()

        this.logFile = task.workDir.resolve(TaskRun.CMD_LOG)
        this.scriptFile = task.workDir.resolve(TaskRun.CMD_SCRIPT)
        this.inputFile =  task.workDir.resolve(TaskRun.CMD_INFILE)
        this.outputFile = task.workDir.resolve(TaskRun.CMD_OUTFILE)
        this.errorFile = task.workDir.resolve(TaskRun.CMD_ERRFILE)
        this.exitFile = task.workDir.resolve(TaskRun.CMD_EXIT)
        this.wrapperFile = task.workDir.resolve(TaskRun.CMD_RUN)
        this.stubFile = task.workDir.resolve(TaskRun.CMD_STUB)
        this.traceFile = task.workDir.resolve(TaskRun.CMD_TRACE)
    }


    /**
     * @return An instance of {@link AwsOptions} holding Batch specific settings
     */
    @Memoized
    protected AwsOptions getAwsOptions() {

        final name = executor.name

        new AwsOptions(
                cliPath: executor.session.getExecConfigProp(name,'awscli',null) as String,
                storageClass: executor.session.config.navigate('aws.client.uploadStorageClass') as String,
                storageEncryption: executor.session.config.navigate('aws.client.storageEncryption') as String,
                remoteBinDir: executor.remoteBinDir as String
            )

    }

    /**
     * {@inheritDoc}
     */
    @Override
    boolean checkIfRunning() {
        if( !jobId )
            return false
        final resp = client.describeJobs(new DescribeJobsRequest().withJobs(jobId))
        if( !resp.getJobs() ) {
            log.debug "[AWS BATCH] cannot retried running status for job=$jobId"
            return false
        }
        final job = resp.getJobs().get(0)
        return job.status in ['RUNNING', 'SUCCEEDED', 'FAILED']
    }

    /**
     * {@inheritDoc}
     */
    @Override
    boolean checkIfCompleted() {
        assert jobId
        final resp = client.describeJobs(new DescribeJobsRequest().withJobs(jobId))
        if( !resp.getJobs() ) {
            log.debug "[AWS BATCH] cannot retried completed status for job=$jobId"
            return false
        }
        final job = resp.getJobs().get(0)
        final done = job.status in ['SUCCEEDED', 'FAILED']
        if( done ) {
            // finalize the task
            task.exitStatus = readExitFile()
            task.stdout = outputFile
            task.stderr = errorFile
            status = TaskStatus.COMPLETED
            return true
        }
        return false
    }

    private int readExitFile() {
        try {
            exitFile.text as Integer
        }
        catch( Exception e ) {
            log.debug "[AWS BATCH] Cannot read exitstatus for task: `$task.name`", e
            return Integer.MAX_VALUE
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    void kill() {
        assert jobId
        log.debug "[AWS BATCH] killing job=$jobId"
        final req = new CancelJobRequest().withJobId(jobId).withReason('Job killed by NF')
        final resp = client.cancelJob(req)
        log.debug "[AWS BATCH] killing job=$jobId; response=$resp"
    }

    /**
     * {@inheritDoc}
     */
    @Override
    void submit() {
        /*
         * create task wrapper
         */
        final launcher = new AwsBatchScriptLauncher(bean,getAwsOptions())
        launcher.scratch = true
        launcher.build()

        final req = newSubmitRequest(task)
        log.trace "[AWS BATCH] now job request > $req"

        /*
         * submit the task execution
         */
        def resp = client.submitJob(req)
        this.jobId = resp.jobId
        this.status = TaskStatus.SUBMITTED
        log.debug "[AWS BATCH] submitted > job=$jobId; work-dir=${task.getWorkDirStr()}"
    }

    /**
     * Retrieve the queue name to use to submit the task execution
     *
     * @param task The {@link TaskRun} object to be executed
     * @return The Batch queue name defined by this job execution
     */
    protected String getJobQueue(TaskRun task) {
        final queue = task.config.queue?.toString()
        if( !queue )
            throw new ProcessNotRecoverableException("Missing AWS Batch job queue -- provide it by using the process `queue` directive")

        return queue
    }

    /**
     * Get the Batch job definition name used to run the specified task
     *
     * @param task The {@link TaskRun} object to be executed
     * @return The Batch job definition name defined by this job execution
     */
    protected String getJobDefinition(TaskRun task) {
        final container = task.getContainer()
        if( !container )
            throw new ProcessNotRecoverableException("Invalid AWS Batch job definition -- provide a Docker image name or a Batch job definition name")

        if( container.startsWith('job-definition://')) {
            return container.substring(17)
        }

        resolveJobDefinition(container)
    }

    /**
     * Maps a docker container image to a Batch job definition name
     *
     * @param container The Docker container image name which need to be used to run the job
     * @return The Batch Job Definition name associated with the specified container
     */
    protected String resolveJobDefinition(String container) {
        if( jobDefinitions.containsKey(container) )
            return jobDefinitions[container]

        // this could race a race condition on the creation of a job definition
        // with another NF instance using the same container name
        synchronized(jobDefinitions) {
            if( jobDefinitions.containsKey(container) )
                return jobDefinitions[container]

            def req = makeJobDefRequest(container)
            def name = findJobDef(req.jobDefinitionName, req.parameters?.nxf_jobid)
            if( name ) {
                log.debug "[AWS BATCH] Found job definition name=$name; container=$container"
            }
            else {
                name = createJobDef(req)
                log.debug "[AWS BATCH] Created job definition name=$name; container=$container"
            }

            jobDefinitions[container] = name
            return name
        }
    }

    /**
     * Create a Batch job definition request object for the specified Docker image
     *
     * @param image The Docker container image for which is required to create a Batch job definition
     * @return An instance of {@link RegisterJobDefinitionRequest} for the specified Docker image
     */
    protected RegisterJobDefinitionRequest makeJobDefRequest(String image) {
        final name = normalizeJobDefinitionName(image)
        final result = new RegisterJobDefinitionRequest()
        result.setJobDefinitionName(name)
        result.setType(JobDefinitionType.Container)

        // container definition
        def container = new ContainerProperties()
                .withImage(image)
        // note the actual command, memory and cpus are overridden when the job is executed
                .withCommand('true')
                .withMemory(1024)
                .withVcpus(1)
        def awscli = getAwsOptions().cliPath
        if( awscli ) {
            def mountName = 'aws-cli'
            def path = Paths.get(awscli).parent.parent.toString()
            def mount = new MountPoint()
                    .withSourceVolume(mountName)
                    .withContainerPath(path)
                    .withReadOnly(true)
            container.setMountPoints([mount])

            def vol = new Volume()
                    .withName(mountName)
                    .withHost(new Host()
                    .withSourcePath(path))
            container.setVolumes([vol])
        }
        result.setContainerProperties(container)

        // create a job marker uuid
        def uuid = CacheHelper.hasher([name, image, awscli]).hash().toString()
        result.setParameters([nxf_jobid:uuid])

        return result
    }

    /**
     * Look for a Batch job definition in ACTIVE status for the given name and NF job definition ID
     *
     * @param name The Batch job definition name
     * @param jobId A unique job definition ID generated by NF
     * @return The fully qualified Batch job definition name eg {@code my-job-definition:3}
     */
    protected String findJobDef(String name, String jobId) {
        log.trace "[AWS BATCH] checking job definition with name=$name; jobid=$jobId"
        final req = new DescribeJobDefinitionsRequest().withJobDefinitionName(name)
        final res = client.describeJobDefinitions(req)
        final jobs = res.getJobDefinitions()
        if( jobs.size()==0 )
            return null

        def job = jobs.find { JobDefinition it -> it.status == 'ACTIVE' && it.parameters?.nxf_jobid == jobId  }
        return job ? "$name:$job.revision" : null
    }

    /**
     * Create (aka register) a new Batch job definition
     *
     * @param req A {@link RegisterJobDefinitionRequest} representing the Batch jib definition to create
     * @return The fully qualified Batch job definition name eg {@code my-job-definition:3}
     */
    protected String createJobDef(RegisterJobDefinitionRequest req) {
        final res = client.registerJobDefinition(req)
        return "${res.jobDefinitionName}:$res.revision"
    }

    /**
     * Make a name string complaint with the Batch job definition format
     *
     * @param name A job name
     * @return A given name formatted to be used as Job definition name
     */
    protected String normalizeJobDefinitionName(String name) {
        if( !name ) return null
        def result = name.replaceAll(/[^a-zA-Z0-9\-_]+/,'-')
        return "nf-" + result
    }

    /**
     * Create a new Batch job request for the given NF {@link TaskRun}
     *
     * @param task A {@link TaskRun} to be executed as Batch job
     * @return A {@link SubmitJobRequest} instance representing the Batch job to submit
     */
    protected SubmitJobRequest newSubmitRequest(TaskRun task) {

        // the cmd list to launch it
        def cli = awsOptions.cliPath ?: 'aws'
        def cmd = ['bash','-c', "$cli s3 cp s3:/${getWrapperFile()} - | bash 2>&1".toString() ] as List<String>

        /*
         * create the request object
         */
        def result = new SubmitJobRequest()
        result.setJobName(normalizeJobName(task.name))
        result.setJobQueue(getJobQueue(task))
        result.setJobDefinition(getJobDefinition(task))

        // set the actual command
        def container = new ContainerOverrides()
        container.command = cmd
        // set the task memory
        if( task.config.getMemory() )
            container.memory = (int)task.config.getMemory().toMega()
        // set the task cpus
        if( task.config.getCpus() > 1 )
            container.vcpus = task.config.getCpus()
        // set the environment
        def vars = getEnvironmentVars()
        if( vars )
            container.setEnvironment(vars)

        result.setContainerOverrides(container)

        return result
    }

    /**
     * @return The list of environment variables to be defined in the Batch job execution context
     */
    protected List<KeyValuePair> getEnvironmentVars() {
        def vars = []
        if( this.environment?.containsKey('NXF_DEBUG') )
            vars << new KeyValuePair().withName('NXF_DEBUG').withValue(this.environment['NXF_DEBUG'])
        return vars
    }

    /**
     * @return The launcher script file {@link Path}
     */
    protected Path getWrapperFile() {
        wrapperFile
    }

    /**
     * Remove invalid characters from a job name string
     *
     * @param name A job name containing possible invalid character
     * @return A job name without invalid characters
     */
    protected String normalizeJobName(String name) {
        name.replaceAll(' ','_').replaceAll(/[^a-zA-Z0-9_]/,'')
    }
}

/**
 * Implements BASH launcher script for AWS Batch jobs
 */
class AwsBatchScriptLauncher extends BashWrapperBuilder {

    AwsBatchScriptLauncher( TaskBean bean, AwsOptions opts ) {
        super(bean, new AwsBatchFileCopyStrategy(bean,opts))
    }

    protected void makeEnvironmentFile(Path file) {
        /**
         * Do nothing here - Environment is managed is built by {@link AwsBatchFileCopyStrategy}
         */
    }
}

/**
 * Defines the script operation to handle file when running in the Cirrus cluster
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class AwsBatchFileCopyStrategy extends SimpleFileCopyStrategy {

    private Path workDir

    private AwsOptions opts

    private Map<String,String> environment

    AwsBatchFileCopyStrategy( TaskBean task, AwsOptions opts ) {
        this.workDir = task.workDir
        this.targetDir = task.getTargetDir()
        this.outputFiles = task.getOutputFiles()
        this.inputFiles = task.getInputFiles()
        // include task script
        this.inputFiles[TaskRun.CMD_SCRIPT] = workDir.resolve(TaskRun.CMD_SCRIPT)
        // include task stdin file
        if( task.input != null ) {
            inputFiles[TaskRun.CMD_INFILE] = workDir.resolve(TaskRun.CMD_INFILE)
        }
        this.opts = opts
        this.environment = task.environment
    }

    /**
     * @return A script snippet that download from S3 the task scripts:
     * {@code .command.env}, {@code .command.sh}, {@code .command.in},
     * etc.
     */
    String getBeforeStartScript() {

        final cli = opts.cliPath ?: 'aws'
        final storage = opts.storageClass ?: 'STANDARD'
		final encryption = opts.storageEncryption ? "--sse $opts.storageEncryption " : ''

        def uploader = """
        # aws helper
        nxf_s3_upload() {
            local pattern=\$1
            local s3path=\$2
            for name in \$pattern;do
              if [[ -d "\$name" ]]; then
                $cli s3 cp \$name \$s3path/\$name --quiet --recursive $encryption--storage-class $storage
              else
                $cli s3 cp \$name \$s3path/\$name --quiet $encryption--storage-class $storage
              fi
          done
        }

        """.stripIndent()

        def env = getEnvScript()
        return env ? uploader + env : uploader
    }

    protected String getEnvScript() {
        final result = new StringBuilder()
        final copy = environment ? new HashMap<String,String>(environment) : Collections.<String,String>emptyMap()
        final path = copy.containsKey('PATH')
        // remove any external PATH
        if( path )
            copy.remove('PATH')
        // when a remote bin directory is provide managed it properly
        if( opts.remoteBinDir ) {
            result << "${opts.cliPath?:'aws'} s3 cp --recursive --quiet s3:/${opts.remoteBinDir} \$PWD/nextflow-bin\n"
            result << "chmod +x \$PWD/nextflow-bin/*\n"
            result << "export PATH=\$PWD/nextflow-bin:\$PATH\n"
        }
        result << TaskProcessor.bashEnvironmentScript(copy)

        return result.size() ? "# process environment\n" + result.toString() + '\n' : null
    }

    /**
     * {@inheritDoc}
     */
    @Override
    String stageInputFile( Path path, String targetName ) {
        final cli = opts.cliPath ?: 'aws'
        def op = "$cli s3 cp --quiet "
        if( path.isDirectory() ) {
            op += "--recursive "
        }
        op += "s3:/${path} ${targetName}"
        return op
    }

    /**
     * {@inheritDoc}
     */
    @Override
    String getUnstageOutputFilesScript() {

        // collect all the expected names (pattern) for files to be un-staged
        def result = []
        def normalized = normalizeGlobStarPaths(outputFiles)

        // create a bash script that will copy the out file to the working directory
        log.trace "[AWS BATCH] Unstaging file path: $normalized"
        if( normalized ) {
            result << ""
            normalized.each {
                result << "nxf_s3_upload '$it' s3:/${targetDir} || true" // <-- add true to avoid it stops on errors
            }
        }

        return result.join(separatorChar)
    }

    /**
     * {@inheritDoc}
     */
    @Override
    String touchFile( Path file ) {
        "touch ${file.name} && nxf_s3_upload ${file.name} s3:/${file.getParent()} && rm ${file.name}"
    }

    /**
     * {@inheritDoc}
     */
    @Override
    String fileStr( Path path ) {
        Escape.path(path.getFileName())
    }

    /**
     * {@inheritDoc}
     */
    @Override
    String copyFile( String name, Path target ) {
        "nxf_s3_upload ${Escape.path(name)} s3:/${Escape.path(target.getParent())}"
    }

    /**
     * {@inheritDoc}
     */
    String exitFile( Path path ) {
        "${path.name} && nxf_s3_upload ${Escape.path(path.name)} s3:/${Escape.path(path.getParent())} || true"
    }

    String pipeInputFile( Path path ) {
        " <(${opts.cliPath ?: 'aws'} s3 cp --quiet s3:/$path -)"
    }

}

/**
 * Helper class wrapping AWS config options required for Batch job executions
 */
@Slf4j
@ToString
@EqualsAndHashCode
@CompileStatic
class AwsOptions {

    String cliPath

    String storageClass

    String storageEncryption

    String remoteBinDir

    void setStorageClass(String value) {
        if( value in [null, 'STANDARD', 'REDUCED_REDUNDANCY'])
            this.storageClass = value
        else
            log.warn "Unsupported AWS storage-class: $value"
    }

    void setStorageEncryption(String value) {
        if( value in [null,'AES256'] )
            this.storageEncryption = value
        else
            log.warn "Unsupported AWS storage-encryption: $value"
    }

    void setCliPath(String value) {
        if( !value )
            this.cliPath = null
        else {
            if( !value.startsWith('/') ) throw new ProcessNotRecoverableException("Not a valid aws-cli tools path: $value -- it must be an absolute path")
            if( !value.endsWith('/bin/aws')) throw new ProcessNotRecoverableException("Not a valid aws-cli tools path: $value -- it must end with the `/bin/aws` suffix")
            this.cliPath = value
        }
    }

}