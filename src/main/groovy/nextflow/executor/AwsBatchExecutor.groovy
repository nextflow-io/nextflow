package nextflow.executor
import java.nio.file.Path

import com.amazonaws.auth.BasicAWSCredentials
import com.amazonaws.regions.RegionUtils
import com.amazonaws.services.batch.AWSBatchClient
import com.amazonaws.services.batch.model.CancelJobRequest
import com.amazonaws.services.batch.model.ContainerOverrides
import com.amazonaws.services.batch.model.DescribeJobsRequest
import com.amazonaws.services.batch.model.KeyValuePair
import com.amazonaws.services.batch.model.SubmitJobRequest
import com.upplication.s3fs.S3Path
import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.Global
import nextflow.exception.AbortOperationException
import nextflow.exception.ProcessNotRecoverableException
import nextflow.processor.TaskBean
import nextflow.processor.TaskHandler
import nextflow.processor.TaskMonitor
import nextflow.processor.TaskPollingMonitor
import nextflow.processor.TaskRun
import nextflow.processor.TaskStatus
import nextflow.util.Duration
import nextflow.util.Escape
/**
 * Experimental AWS Batch executor
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class AwsBatchExecutor extends Executor {

    @PackageScope
    static AWSBatchClient client

    @Override
    void register() {
        super.register()

        /*
         * make sure the work dir is a S3 bucket
         */
        if( !(session.workDir instanceof S3Path) ) {
            session.abort()
            throw new AbortOperationException("When using `aws-batch` executor a S3 bucket must be provided as working directory -- Add the option `-w s3://<your-bucket/path>` to your run command line")
        }

        /*
         * retrieve config and credentials and create AWS client
         */
        def keys = Global.getAwsCredentials()
        if( !keys )
            throw new IllegalStateException("Missing AWS credentials")

        // TODO improve to remove deprecated API
        def credentials = new BasicAWSCredentials(keys[0], keys[1])
        this.client = new AWSBatchClient(credentials)

        def str = Global.getAwsRegion()
        if( str ) {
            def region = RegionUtils.getRegion(Global.getAwsRegion());
            if( region == null )
                throw new IllegalArgumentException("Not a valid S3 region name: ${str}");
            client.setRegion(region);
        }

    }

    @Override
    protected TaskMonitor createTaskMonitor() {
        TaskPollingMonitor.create(session, name, 50, Duration.of('5 sec'))
    }

    @Override
    TaskHandler createTaskHandler(TaskRun task) {
        assert task
        assert task.workDir

        log.debug "Launching process > ${task.name} -- work folder: ${task.workDir}"
        new AwsBatchTaskHandler(task, this)
    }
}


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


    /** only for testing purpose -- do not use */
    protected AwsBatchTaskHandler() {}


    AwsBatchTaskHandler(TaskRun task, AwsBatchExecutor executor) {
        super(task)
        this.bean = new TaskBean(task)
        this.executor = executor
        this.client = executor.client

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
            log.debug "Cannot read exitstatus for task: `$task.name`", e
            return Integer.MAX_VALUE
        }
    }

    @Override
    void kill() {
        assert jobId
        log.debug "[AWS BATCH] killing job=$jobId"
        final req = new CancelJobRequest().withJobId(jobId).withReason('Job killed by NF')
        final resp = client.cancelJob(req)
        log.debug "[AWS BATCH] killing job=$jobId; response=$resp"
    }

    @Override
    void submit() {

        /*
         * create task wrapper
         */
        def launcher = new AwsBatchScriptLauncher(bean)
        launcher.scratch = true
        launcher.build()

        final req = newSubmitRequest(task)
        log.debug "[AWS BATCH] now job request > $req"

        /*
         * submit the task execution
         */
        def resp = client.submitJob(req)
        this.jobId = resp.jobId
        this.status = TaskStatus.SUBMITTED
        log.debug "[AWS BATCH] submitted > job=$jobId; work-dir=${task.getWorkDirStr()}"
    }

    protected List<String> getJobQueueAndDefinition(TaskRun task) {
        def result = task.config.queue?.toString()
        if( !result )
            throw new ProcessNotRecoverableException("Missing AWS Batch job queue -- provide it by using the process `queue` directive")

        if( !result.contains('/') )
            throw new ProcessNotRecoverableException("Invalid AWS Batch job queue definition -- it must be in the form: jobQueue/jobDefinition")

        result.tokenize('/')
    }

    protected SubmitJobRequest newSubmitRequest(TaskRun task) {

        def awsCli = Global.GetAwsCliPath() 
        // the cmd list to launch it
        def cmd = ['bash','-c', "$awsCli s3 cp s3:/${wrapperFile} - | bash 2>&1 | tee $TaskRun.CMD_LOG".toString() ] as List<String>

        def jobQueueAndDefinition = getJobQueueAndDefinition(task)
        /*
         * create the request object
         */
        def result = new SubmitJobRequest()
        result.setJobName(normalizeName(task.name))
        result.setJobQueue(jobQueueAndDefinition[0])
        result.setJobDefinition(jobQueueAndDefinition[1])

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
        def vars = []
        if( bean.environment )
            bean.environment.each { k,v -> vars << new KeyValuePair().withName(k).withValue(v) }
        if( System.getenv('NXF_DEBUG') )
            vars << new KeyValuePair().withName('NXF_DEBUG').withValue(System.getenv('NXF_DEBUG'))
        container.setEnvironment(vars)

        result.setContainerOverrides(container)

        return result
    }

    protected String normalizeName(String str) {
        str.replaceAll(' ','_').replaceAll(/[^a-zA-Z0-9_]/,'')
    }
}


class AwsBatchScriptLauncher extends BashWrapperBuilder {

    AwsBatchScriptLauncher( TaskBean bean ) {
        super(bean, new AwsBatchFileCopyStrategy(bean))
    }

    protected boolean containerInit() { false }

    protected void makeEnvironmentFile(Path file) { /* do nothing */  }
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

    AwsBatchFileCopyStrategy( TaskBean task ) {
        workDir = task.workDir
        targetDir = task.getTargetDir()
        outputFiles = task.getOutputFiles()
        inputFiles = task.getInputFiles()
        // include task script
        inputFiles[TaskRun.CMD_SCRIPT] = workDir.resolve(TaskRun.CMD_SCRIPT)
        // include task stdin file
        if( task.input != null ) {
            inputFiles[TaskRun.CMD_INFILE] = workDir.resolve(TaskRun.CMD_INFILE)
        }
    }


    /**
     * @return A script snippet that download from S3 the task scripts:
     * {@code .command.env}, {@code .command.sh}, {@code .command.in},
     * etc.
     */
    String getBeforeStartScript() {
	
        def awsCli = Global.GetAwsCliPath() 
		def S3StorageClass = Global.AwsGetStorageClass()
		def S3Encryption = Global.AwsGetStorageEncryption()
		def S3EncryptionCommand = ""
		if(S3Encryption == "AES256") {
			S3EncryptionCommand += "--sse "+S3Encryption
		}

        """
        # aws helper
        nxf_s3_upload() {
            local pattern=\$1
            local s3path=\$2
            for name in \$pattern;do
              if [[ -d "\$name" ]]; then
                echo "$awsCli s3 cp \$name \$s3path/\$name --quiet --storage-class $S3StorageClass --recursive $S3EncryptionCommand"
                $awsCli s3 cp \$name \$s3path/\$name --quiet --storage-class $S3StorageClass --recursive $S3EncryptionCommand
              else
                echo "$awsCli s3 cp \$name \$s3path/\$name --quiet --storage-class $S3StorageClass $S3EncryptionCommand"
                $awsCli s3 cp \$name \$s3path/\$name --quiet --storage-class $S3StorageClass $S3EncryptionCommand
              fi
          done
        }

        """.stripIndent()

    }

    /**
     * {@inheritDoc}
     */
    @Override
    String stageInputFile( Path path, String targetName ) {
        def op = "aws s3 cp --quiet "
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
        log.debug "Unstaging file path: $normalized"
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
}
