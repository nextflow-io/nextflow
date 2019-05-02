/*
 * Copyright 2019, Genome Research Limited
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

package nextflow.executor

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import java.nio.file.Path
import java.nio.file.Paths

import java.security.GeneralSecurityException;
import java.security.KeyStore;
import java.security.SecureRandom;
import java.security.cert.Certificate;
import java.security.cert.CertificateException;
import java.security.cert.CertificateFactory;
import java.security.cert.X509Certificate;
import javax.net.ssl.HostnameVerifier;
import javax.net.ssl.KeyManager;
import javax.net.ssl.SSLContext;
import javax.net.ssl.SSLSession;
import javax.net.ssl.TrustManager;
import javax.net.ssl.TrustManagerFactory;
import javax.net.ssl.X509TrustManager;
import javax.net.ssl.HttpsURLConnection

import nextflow.exception.AbortOperationException
import nextflow.processor.TaskHandler
import nextflow.processor.TaskMonitor
import nextflow.processor.WrMonitor
import nextflow.processor.TaskRun
import nextflow.processor.TaskBean
import nextflow.processor.TaskStatus
import nextflow.processor.BatchContext
import nextflow.processor.BatchHandler
import nextflow.util.Escape
import groovy.json.JsonSlurper
import static groovy.json.JsonOutput.toJson
/**
 * Executor that schedules jobs using wr as a backend, avoiding storage of
 * state on disk.
 *
 * See https://github.com/VertebrateResequencing/wr
 *
 * @author Sendu Bala <sb10@sanger.ac.uk>
 * Based on TesExecutor by Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class WrExecutor extends Executor {

    static private String token
    static private String defaultManagerDir
    static private String endpoint
    static private String cacertPath
    static private WrRestApi client

    @Override
    void register() {
        if( session.binDir && !session.binDir.empty() ) {
            session.abort()
            throw new AbortOperationException("ERROR: wr executor does not allow the use of custom scripts in the `bin` folder")
        }

        super.register()
    }

    protected String getDisplayName() {
        return "$name [$endpoint]"
    }

    WrRestApi getClient() {
        if (!client) {
            defaultManagerDir = Paths.get(System.getProperty('user.home'), ".wr_production")
            endpoint = getEndPoint()
            client = new WrRestApi(endpoint, getToken(), getCacertPath())
        }
        client
    }

    protected String getEndPoint() {
        // default port that wr listens on is 1021 + (uid * 4) + 1
        // *** note that this will probably only work on linux/mac os, but wr
        // probably only works fully on those as well...
        int uid = ["id", "-u"].execute().text.trim() as Integer
        int port = 1021 + (uid * 4) + 1

        def result = session.getConfigAttribute('executor.endpoint', "https://localhost:$port")
        log.debug "[wr] endpoint=$result"
        return result
    }

    protected String getToken() {
        String path = session.getConfigAttribute('executor.tokenpath', Paths.get(defaultManagerDir, "client.token"))
        log.debug "[wr] tokenpath=$path"
        String result = new File(path).text
        return result
    }

    protected String getCacertPath() {
        String path = session.getConfigAttribute('executor.cacertpath', Paths.get(defaultManagerDir, "ca.pem"))
        log.debug "[wr] cacertpath=$path"
        return path
    }

    /**
     * @return {@code false} whenever the containerization is managed by the executor itself
     */
    boolean isContainerNative() {
        return false
    }

    /**
     * Create a queue holder for this executor
     *
     * @return
     */
    TaskMonitor createTaskMonitor() {
        return WrMonitor.create(session, getClient())
    }

    /*
     * Prepare and launch the task in the underlying execution platform
     */
    @Override
    TaskHandler createTaskHandler(TaskRun task) {
        assert task
        assert task.workDir
        log.debug "[wr] Launching process > ${task.name} -- work folder: ${task.workDir}"
        new WrTaskHandler(task, this)
    }
    
}


/**
 * Handles a job execution using wr as a backend, without needing files on disk
 * for state.
 *
 * @author Sendu Bala <sb10@sanger.ac.uk>
 * Based on TesTaskHandler by Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class WrTaskHandler extends TaskHandler implements BatchHandler<String,Map> {

    public static final List<String> COMPLETE_STATUSES = ['complete', 'buried', 'deleted']
    public static final List<String> STARTED_STATUSES = ['delayed', 'reserved', 'running'] + COMPLETE_STATUSES

    final WrExecutor executor
    private WrRestApi client
    private final Path wrapperFile
    private final Path outputFile
    private final Path errorFile
    private final Path logFile
    private final Path scriptFile
    private final Path inputFile
    private final Path traceFile
    private String jobId

    /**
     * Batch context shared between multiple task handlers
     */
    private BatchContext<String,Map> context

    WrTaskHandler(TaskRun task, WrExecutor executor) {
        super(task)
        this.executor = executor
        this.client = executor.getClient()

        this.logFile = task.workDir.resolve(TaskRun.CMD_LOG)
        this.scriptFile = task.workDir.resolve(TaskRun.CMD_SCRIPT)
        this.inputFile =  task.workDir.resolve(TaskRun.CMD_INFILE)
        this.outputFile = task.workDir.resolve(TaskRun.CMD_OUTFILE)
        this.errorFile = task.workDir.resolve(TaskRun.CMD_ERRFILE)
        this.wrapperFile = task.workDir.resolve(TaskRun.CMD_RUN)
        this.traceFile = task.workDir.resolve(TaskRun.CMD_TRACE)
    }

    void batch( BatchContext<String,Map> context ) {
        if( jobId ) {
            context.collect(jobId)
            this.context = context
        }
    }

    private String jobIdsToString(Collection batchIds) {
        final MAX=10
        final sz = batchIds.size()
        batchIds.size()<=MAX ? batchIds.join(',').toString() : batchIds.take(MAX).join(',').toString() + ", ... other ${sz-MAX} omitted"
    }

    /**
     * Retrieve Batch job status information
     *
     * @param jobId The Batch job ID
     * @return The associated job details in a Map or {@code null} if no information is found
     */
    protected Map getJob(String jobId) {
        Collection batchIds
        if( context ) {
            // check if this response is cached in the batch collector
            if( context.contains(jobId) ) {
                log.trace "[wr] hit cache when getting job=$jobId"
                return context.get(jobId)
            }
            log.trace "[wr] missed cache when getting job=$jobId"
            batchIds = context.getBatchFor(jobId, 1000)
        }
        else {
            batchIds = [jobId]
        }

        // retrieve the status for the specified batch of jobs
        final String ids = jobIdsToString(batchIds)
        log.trace "[wr] getting jobs=${jobIdsToString(batchIds)}"
        List<Map> jobs = client.status(batchIds.join(',').toString())
        if( !jobs ) {
            log.debug "[wr] cannot retrieve status for job=$jobId"
            return null
        }

        Map result=null
        jobs.each {
            String id = it."Key" as String
            // cache the response in the batch collector
            context?.put( id, it )
            // return the job detail for the specified job
            if( id == jobId )
                result = it
        }
        if( !result ) {
            log.debug "[wr] cannot find status for job=$jobId"
        }

        return result
    }

    @Override
    boolean checkIfRunning() {
        if( !jobId || !isSubmitted() )
            return false
        final job = getJob(jobId)
        final result = job?.State in STARTED_STATUSES
        if( result ) {
            log.trace "[wr] Task started > $task.name"
            this.status = TaskStatus.RUNNING
        }
        return result
    }

    @Override
    boolean checkIfCompleted() {
        if( !isRunning() )
            return false
        final job = getJob(jobId)
        final done = job?.State in COMPLETE_STATUSES
        if( done ) {
            // finalize the task
            log.trace "[wr] Task completed > $task.name"
            if (job.Exited) {
                task.exitStatus = job.Exitcode as Integer
            } else {
                task.exitStatus = Integer.MAX_VALUE
            }
            task.stdout = outputFile
            task.stderr = errorFile
            status = TaskStatus.COMPLETED
        }
        return done
    }

    @Override
    void kill() {
        if( jobId )
            client.cancel(jobId)
        else
            log.trace "[wr] Invalid kill request -- missing jobId"
    }

    @Override
    void submit() {
        // this is not actually called by anything since WrMonitor does batched
        // submits direcly to the client, but this is implemented anyway
        List<Map> jobs = client.add([submitArgs()])
        jobId = jobs[0]."Key" as String
        // log.debug("[wr] submitted job $jobId")
        status = TaskStatus.SUBMITTED
    }

    List submitArgs() {
        // create task wrapper
        final bash = new WrBashBuilder(task)
        bash.build()

        WrFileCopyStrategy copyStrategy = bash.copyStrategy as WrFileCopyStrategy
        String wrapperPath = copyStrategy.wrWorkPath(wrapperFile)

        return ["/bin/bash $wrapperPath", task, copyStrategy]
    }

    void submitted(String id) {
        jobId = id
        status = TaskStatus.SUBMITTED
    }

}

/**
 * WrRestApi lets you interact with wr's REST API.
 *
 * See https://github.com/VertebrateResequencing/wr/wiki/REST-API
 *
 * @author Sendu Bala <sb10@sanger.ac.uk>
 */
@Slf4j
@CompileStatic
class WrRestApi {
    
    static private String endpoint
    static private String token

    static private final String rest = "/rest/v1"
    static private final String jobs = "/jobs/"

    WrRestApi(String endpoint, String token, String cacertPath) {
        this.endpoint = endpoint
        this.token = token

        // have Groovy trust cacertPath
        def nullHostnameVerifier = [
            verify: { hostname, session -> true }
        ]
        HttpsURLConnection.setDefaultHostnameVerifier(nullHostnameVerifier as HostnameVerifier)

        def is = new File(cacertPath).newInputStream()

        TrustManager[] trustManagers = null;
        char[] password = null; // any password will work.
        KeyStore caKeyStore = KeyStore.getInstance(KeyStore.getDefaultType());
        caKeyStore.load(null, password);
        CertificateFactory certificateFactory = CertificateFactory.getInstance("X.509");
        Collection<? extends Certificate> certificates = certificateFactory.generateCertificates(is);
        if (certificates.isEmpty()) {
            throw new IllegalArgumentException("expected non-empty set of trusted certificates");
        }
        int index = 0;
        certificates.each {
            String certificateAlias = "ca" + Integer.toString(index++);
            caKeyStore.setCertificateEntry(certificateAlias, it);
        }
        TrustManagerFactory trustManagerFactory = TrustManagerFactory.getInstance(TrustManagerFactory.getDefaultAlgorithm());
        trustManagerFactory.init(caKeyStore);
        trustManagers = trustManagerFactory.getTrustManagers();

        SSLContext sslContext = SSLContext.getInstance("SSL");
        sslContext.init(null, trustManagers, new SecureRandom());
        HttpsURLConnection.setDefaultSSLSocketFactory(sslContext.getSocketFactory());

        is.close()
    }

    /**
     * Adds new commands to the job queue. Takes a List of Lists, where the
     * sub lists are [String cmd, TaskRun task, WrFileCopyStrategy copyStrategy]
     * Returns the jobs.
     */
    List<Map> add(List<List> argList) {
        // curl --cacert ~/.wr_production/ca.pem -H "Content-Type: application/json" -H "Authorization: Bearer [token]" -X POST -d '[{"cmd":"sleep 5 && echo mymsg && false","memory":"5M","cpus":1},{"cmd":"sleep 5","cpus":1}]' 'https://localhost:11302/rest/v1/jobs/?rep_grp=myid&cpus=2&memory=3G&time=5s'

        List<Map> jsonArgs = []

        argList.each {
            String cmd = it[0]
            TaskRun task = it[1]
            WrFileCopyStrategy copyStrategy = it[2]

            String cwd = ""
            boolean cwdMatters = true
            switch (copyStrategy.getPathScheme(task.workDir)) {
                case "file":
                    cwd = task.getWorkDirStr()
                    break
                case "s3":
                    cwd = "/tmp"
                    cwdMatters = false
                    break
                default:
                    throw new IllegalArgumentException("Unsupported scheme for work dir: ${task.getWorkDirStr()}")
            }

            String grp = "[nextflow] $task.processor.name"

            // calculate mounts option
            def reads = copyStrategy.inputBuckets()
            List<Map> m = []
            for (path in reads) {
                m << ["Mount":copyStrategy.getInputMountPath() + "/" + path, "Targets":[["Path":path]]]
            }
            def write = copyStrategy.outputBucket()
            if (write) {
                m << ["Mount":copyStrategy.getOutputMountPath(),"Targets":[["Path":write,"Write":true]]]
            }

            Integer cpus = 1
            if ( task.config.cpus && task.config.cpus > 0 ) {
                cpus = task.config.cpus
            }

            Integer override = 0

            String mem = ""
            if( task.config.getMemory() ) {
                mem = String.valueOf(task.config.getMemory().toUnit('MB')) + "M"
                override = 2
            }

            String t = ""
            if( task.config.time ) {
                t = task.config.getTime().format('Hh')
                override = 2
            }

            Long d = 0
            if( task.config.getDisk() ) {
                def disk = task.config.getDisk()
                d = disk.toUnit('GB')
                override = 2
            }

            List<String> limits = []
            if ( task.config.maxForks ) {
                limits << sprintf('%s:%d', task.processor.name, task.config.maxForks)
            }

            // *** what about cloud opts like image and flavor?
            // Turn on docker monitoring if docker container is being used?
            // Does BashWrapperBuilder handle everything to do with env vars?

            Map args = [cmd: cmd, cwd: cwd, cwd_matters: cwdMatters, rep_grp: grp, req_grp: grp, limit_grps: limits, override: override, retries: 0, cpus: cpus, memory: mem, time: t, disk: d, mounts: m]
            // log.debug "[wr] add args: $args"
            jsonArgs << args
        }

        def response = postJson(jobs, jsonArgs)
        // log.debug("made a POST request")
        return parseJobsFromJson(response)
    }

    /**
     * Get the status of a previously add()ed job(s). Returns a List of Map with
     * details of each job, including Key, State, Exited and Exitcode.
     */
    List<Map> status(String ids) {
        // curl --cacert ~/.wr_production/ca.pem -H "Authorization: Bearer [token]" 'https://localhost:11302/rest/v1/jobs/58cef10e7a340c3b7fa09ea304a3cb98'

        // log.debug "working on ids $ids"

        def get = authenticatedConnection("$jobs/$ids")
        get.setRequestMethod("GET")
        def getRC = get.getResponseCode()
        // log.debug("made a GET request")
        if (!getRC.equals(200)) {
            log.debug("[wr] failed to get command details from wr manager")
            return null
        }
        return parseJobsFromJson(get.getInputStream().getText())
    }

    /**
     * Cancel the given job. If it's running it will be killed. If it's lost
     * it will be buried. If it's otherwise incomplete, it will be removed from
     * wr's queue.
     */
    void cancel(String id) {
        // like status, but send DELETE instead of GET. Also state parameter
        // required to be running|lost|deletable

        // *** is there a way to avoid needing to get the current state of the
        // job first? Does our caller only call us for running jobs?
        List jobs = status(id)
        Map job = jobs[0]
        String state
        if (job.State == 'running') {
            state = 'running'
        } else if (job.State == 'lost') {
            state = 'lost'
        } else if (job.State != 'complete') {
            state = 'deletable'
        }

        def delete = authenticatedConnection("$jobs/$id?state=$state")
        delete.setRequestMethod("DELETE")
        def deleteRC = delete.getResponseCode()
        if (!deleteRC.equals(200) || !deleteRC.equals(202)) {
            log.debug("[wr] failed to cancel command held in wr manager")
        }
    }

    private HttpsURLConnection authenticatedConnection(String leaf) {
        String url = "$endpoint$rest$leaf"
        // log.debug("[wr] url: $url")
        def post = new URL(url).openConnection() as HttpsURLConnection
        post.setRequestProperty("Content-Type", "application/json")
        post.setRequestProperty("Authorization", "Bearer $token")
        return post
    }

    private HttpsURLConnection postConnection(String leaf) {
        def post = authenticatedConnection(leaf)
        post.setRequestMethod("POST")
        post.setDoOutput(true)
        return post
    }

    private String postJson(String leaf, List<Map> args) {
        def json = toJson(args)
        // log.debug("[wr] posting JSON: $json")

        def post = postConnection(leaf)
        post.getOutputStream().write(json.getBytes("UTF-8"))
        def postRC = post.getResponseCode()
        // log.debug(postRC as String)
        if (!postRC.equals(201)) {
            log.debug("[wr] failed to send new command details to wr manager")
            return null
        }
        return post.getInputStream().getText()
    }

    private List<Map> parseJobsFromJson(String json) {
        List<Map> jobs = new JsonSlurper().parseText(json) as List<Map>
        // log.debug "got ${jobs.size()} jobs"
        return jobs
    }

    private String parseIdFromJson(String json) {
        return extractFirstJob(parseJobsFromJson(json)).Key
    }

    private Map extractFirstJob(List<Map> jobs) {
        if (jobs.size() != 1) {
            log.debug("[wr] expected 1 job, got ${jobs.size()}")
        }
        Map job = jobs[0]
        return job
    }

}

/**
 * Bash builder adapter to manage wr specific tasks
 *
 * @author Sendu Bala <sb10@sanger.ac.uk>
 * based on TesBashBuilder by Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class WrBashBuilder extends BashWrapperBuilder {

    WrBashBuilder(TaskRun task) {
        super(new TaskBean(task), new WrFileCopyStrategy(new TaskBean(task)))
    }

    WrBashBuilder(TaskBean bean) {
        super(bean, new WrFileCopyStrategy(bean))
    }

    protected boolean alwaysTryToUnstage() {
        return true
    }

}

/**
 * Defines the script operation to handle files when running via wr.
 *
 * @author Sendu Bala <sb10@sanger.ac.uk>
 * based on AwsBatchFileCopyStrategy by Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class WrFileCopyStrategy extends SimpleFileCopyStrategy {

    private final String outputMountLocation
    static final String inputMountLocation = ".inputs"

    private Map<String,String> environment

    private final Map<String,Path> inputFiles

    WrFileCopyStrategy(TaskBean task) {
        super(task)
        this.environment = task.environment
        this.inputFiles = task.inputFiles
        this.outputMountLocation = ".mnt" + task.workDir.toString()
    }

    /**
     * {@inheritDoc}
     */
    String getBeforeStartScript() {
        return null
    }

    /**
     * {@inheritDoc}
     */
    @Override
    String getEnvScript(Map environment, boolean container) {
        if( container )
            throw new IllegalArgumentException("Parameter `wrapHandler` not supported by ${this.class.simpleName}")

        final result = new StringBuilder()
        final copy = environment ? new HashMap<String,String>(environment) : Collections.<String,String>emptyMap()
        final path = copy.containsKey('PATH')
        // remove any external PATH
        if( path )
            copy.remove('PATH')
        // when a remote bin directory is provided manage it properly
        // *** somehow accept a remoteBinDir opt, and assume it will be mounted
        // somewhere, and add that relative path to PATH?...
        // if( opts.remoteBinDir ) {
        //     result << "${opts.getAwsCli()} s3 cp --recursive --only-show-errors s3:/${opts.remoteBinDir} \$PWD/nextflow-bin\n"
        //     result << "chmod +x \$PWD/nextflow-bin/*\n"
        //     result << "export PATH=\$PWD/nextflow-bin:\$PATH\n"
        // }
        // finally render the environment
        final envSnippet = super.getEnvScript(copy,false)
        if( envSnippet )
            result << envSnippet
        return result.toString()
    }

    /**
     * wrWorkPath returns the same paths that SimpleFileCopyStrategy would end
     * up using, unless path is in S3, in which case it is the file path in the
     * output mount location.
    */
    String wrWorkPath( Path path ) {
        return wrPath(path, outputMountLocation)
    }

    /**
     * wrInputPath returns the same paths that SimpleFileCopyStrategy would end
     * up using, unless path is in S3, in which case it is the file path in the
     * input mount location.
    */
    String wrInputPath( Path path ) {
        return wrPath(path, inputMountLocation)
    }

    /**
     * wrPath returns the same paths that SimpleFileCopyStrategy would end
     * up using, unless path is in S3, in which case it is the file path in the
     * given mount location.
    */
    private String wrPath( Path path, String mountLocation ) {
        if( getPathScheme(path) == 's3' ) {
            if (path == workDir) {
                return "$mountLocation/"
            }
            if (mountLocation == outputMountLocation) {
                String baseName = Escape.path(path.getFileName())
                return "$mountLocation/$baseName"
            }
            String p = Escape.path(path)
            return "$mountLocation$p"
        }
        return Escape.path(path)
    }

    /**
     * getInputMountPath tells you where you should mount the S3 locations
     * returned by inputBuckets(). Each different bucket should be mounted on
     * its own subdirectoy of this path, named after the bucket.
    */
    String getInputMountPath() {
        return inputMountLocation
    }

    /**
     * getOutputMountPath tells you where you should mount the S3 location
     * returned by outputBucket().
    */
    String getOutputMountPath() {
        return outputMountLocation
    }

    /**
     * inputBuckets returns all the different S3 buckets that would need to be
     * mounted at getInputMountPath() for access to the input files.
    */
    List<String> inputBuckets() {
        List<String> result = []
        inputFiles.each{
            if (getPathScheme(it.value) == 's3') {
                String dir = it.value.getRoot().toString().substring(1)
                dir = dir.substring(0, dir.length() - 1)
                if (!result.contains(dir)) {
                    result << dir
                }
            }
        }
        return result
    }

    /**
     * outputBucket returns the S3 bucket directory that would need to be
     * mounted at getOutputMountPath() in order to output files.
    */
    String outputBucket() {
        if (getPathScheme(workDir) == "s3") {
            return workDir.toString().substring(1)
        }
        return ""
    }

    /**
     * {@inheritDoc}
     */
    @Override
    String stageInputFile( Path path, String targetName ) {
        def cmd = ''
        def p = targetName.lastIndexOf('/')
        if( p>0 ) {
            cmd += "mkdir -p ${Escape.path(targetName.substring(0,p))} && "
        }
        cmd += stageInCommand(wrInputPath(path.toAbsolutePath()), targetName, stageinMode)

        return cmd
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected String stageOutCommand( String source, Path targetDir, String mode ) {
        return stageOutCommand(source, wrWorkPath(targetDir), mode)
    }

    /**
     * {@inheritDoc}
     */
    @Override
    String touchFile( Path path ) {
        "touch ${wrWorkPath(path)}"
    }

    /**
     * {@inheritDoc}
     */
    @Override
    String fileStr( Path path ) {
        wrWorkPath(path)
    }

    /**
     * {@inheritDoc}
     */
    @Override
    String copyFile( String name, Path target ) {
        "cp ${Escape.path(name)} ${wrWorkPath(target)}"
    }

    /**
     * {@inheritDoc}
     */
    String exitFile( Path path ) {
        "> ${wrWorkPath(path)}"
    }

    /**
     * {@inheritDoc}
     */
    @Override
    String pipeInputFile( Path path ) {
        " < ${wrInputPath(path)}"
    }
}