/*
 * Copyright 2013-2024, Seqera Labs
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
 *
 */

package io.seqera.tower.plugin


import java.nio.file.Path
import java.time.Instant
import java.time.OffsetDateTime
import java.time.ZoneId
import java.util.concurrent.LinkedBlockingQueue
import java.util.concurrent.TimeUnit

import groovy.json.JsonGenerator
import groovy.json.JsonOutput
import groovy.json.JsonSlurper
import groovy.transform.CompileStatic
import groovy.transform.ToString
import groovy.transform.TupleConstructor
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.exception.AbortOperationException
import nextflow.processor.TaskHandler
import nextflow.processor.TaskId
import nextflow.processor.TaskProcessor
import nextflow.trace.ResourcesAggregator
import nextflow.trace.TraceObserver
import nextflow.trace.TraceRecord
import nextflow.util.Duration
import nextflow.util.LoggerHelper
import nextflow.util.ProcessHelper
import nextflow.util.SimpleHttpClient
import nextflow.util.Threads

/**
 * Send out messages via HTTP to a configured URL on different workflow
 * execution events.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class TowerClient implements TraceObserver {

    static final public String DEF_ENDPOINT_URL = 'https://api.cloud.seqera.io'

    static private final int TASKS_PER_REQUEST = 100

    static private final Duration REQUEST_INTERVAL = Duration.of('1 sec')

    static private final Duration ALIVE_INTERVAL = Duration.of('1 min')

    static private final String TOKEN_PREFIX = '@token:'

    @TupleConstructor
    static class Response {
        final int code
        final String message
        final String cause
        boolean isError() { code < 200 || code >= 300 }
    }

    @ToString(includeNames = true)
    static class ProcessEvent {
        TraceRecord trace
        boolean completed
    }
    

    private Session session

    /**
     * Workflow identifier, will be taken from the Session() object later
     */
    private String runName

    /**
     * Store the sessions unique ID for downstream reference purposes
     */
    private String runId

    /**
     * Simple http client object that will send out messages
     */
    private SimpleHttpClient httpClient

    private JsonGenerator generator

    private String workflowId

    private String watchUrl

    private String endpoint

    private ResourcesAggregator aggregator

    protected Map<String,String> env = System.getenv()

    private LinkedBlockingQueue<ProcessEvent> events = new LinkedBlockingQueue()

    private Thread sender

    private Duration requestInterval = REQUEST_INTERVAL

    private Duration aliveInterval = ALIVE_INTERVAL

    private LinkedHashSet<String> processNames = new LinkedHashSet<>(20)

    private boolean terminated

    private Map<String,Integer> schema = Collections.emptyMap()

    private int maxRetries = 5

    private int backOffDelay

    private int backOffBase

    private boolean towerLaunch

    private String refreshToken

    private String workspaceId

    private TowerReports reports


    /**
     * Constructor that consumes a URL and creates
     * a basic HTTP client.
     * @param endpoint The target address for sending messages to
     */
    TowerClient(Session session, String endpoint) {
        this.session = session
        this.endpoint = checkUrl(endpoint)
        this.schema = loadSchema()
        this.generator = TowerJsonGenerator.create(schema)
        this.reports = new TowerReports(session)
    }

    TowerClient withEnvironment(Map env) {
        this.env = env
        return this
    }

    /**
     * only for testing purpose -- do not use
     */
    protected TowerClient() {
        this.generator = TowerJsonGenerator.create(Collections.EMPTY_MAP)
    }

    boolean enableMetrics() { true }

    String getEndpoint() { endpoint }

    String getWorkflowId() { workflowId }

    boolean getTowerLaunch() { towerLaunch }

    String getRunName() { runName }

    String getRunId() { runId }

    void setAliveInterval(Duration d) {
        this.aliveInterval = d
    }

    void setRequestInterval(Duration d) {
        this.requestInterval = d
    }

    void setMaxRetries( int value ) {
        this.maxRetries = value
    }

    void setBackOffBase( int value ) {
        this.backOffBase = value
    }

    void setBackOffDelay( int value ) {
        this.backOffDelay = value
    }

    void setWorkspaceId( String workspaceId ) {
        this.workspaceId = workspaceId
    }

    String getWorkspaceId() { workspaceId }

    /**
     * Check the URL and create an HttpPost() object. If a invalid i.e. protocol is used,
     * the constructor will raise an exception.
     *
     * The RegEx was taken and adapted from http://urlregex.com
     *
     * @param url String with target URL
     * @return The requested url or the default url, if invalid
     */
    protected String checkUrl(String url){
        // report a warning for legacy endpoint
        if( url.contains('https://api.tower.nf') ) {
            log.warn "The endpoint `https://api.tower.nf` is deprecated - Please use `https://api.cloud.seqera.io` instead"
        }
        if( url =~ "^(https|http)://[-a-zA-Z0-9+&@#/%?=~_|!:,.;]*[-a-zA-Z0-9+&@#/%=~_|]" ) {
            while( url.endsWith('/') )
                url = url[0..-2]
            return url
        }
        throw new IllegalArgumentException("Only http and https are supported -- The given URL was: ${url}")
    }

    protected String getHostUrl(String endpoint) {
        def url = new URL(endpoint)
        return "${url.protocol}://${url.authority}"
    }

    protected String getUrlTraceCreate() {
        def result = this.endpoint + '/trace/create'
        if( workspaceId )
            result += "?workspaceId=$workspaceId"
        return result
    }

    protected String getUrlTraceBegin() {
        def result = "$endpoint/trace/$workflowId/begin"
        if( workspaceId )
            result += "?workspaceId=$workspaceId"
        return result
    }

    protected String getUrlTraceComplete() {
        def result = "$endpoint/trace/$workflowId/complete"
        if( workspaceId )
            result += "?workspaceId=$workspaceId"
        return result
    }

    protected String getUrlTraceHeartbeat() {
        def result = "$endpoint/trace/$workflowId/heartbeat"
        if( workspaceId )
            result += "?workspaceId=$workspaceId"
        return result
    }

    protected String getUrlTraceProgress() {
        def result = "$endpoint/trace/$workflowId/progress"
        if( workspaceId )
            result += "?workspaceId=$workspaceId"
        return result
    }

    /**
     * On workflow start, submit a message with some basic
     * information, like Id, activity and an ISO 8601 formatted
     * timestamp.
     * @param session The current Nextflow session object
     */
    @Override
    void onFlowCreate(Session session) {
        log.debug "Creating Seqera Platform observer -- endpoint=$endpoint; requestInterval=$requestInterval; aliveInterval=$aliveInterval; maxRetries=$maxRetries; backOffBase=$backOffBase; backOffDelay=$backOffDelay"
        
        this.session = session
        this.aggregator = new ResourcesAggregator(session)
        this.runName = session.getRunName()
        this.runId = session.getUniqueId()
        this.httpClient = new SimpleHttpClient()
        // set the auth token
        setAuthToken( httpClient, getAccessToken() )

        // send hello to verify auth
        final req = makeCreateReq(session)
        final resp = sendHttpMessage(urlTraceCreate, req, 'POST')
        if( resp.error ) {
            log.debug """\
                Unexpected HTTP response
                - endpoint    : $urlTraceCreate
                - status code : $resp.code
                - response msg: $resp.cause
                """.stripIndent(true)
            throw new AbortOperationException(resp.message)
        }
        final ret = parseTowerResponse(resp)
        this.workflowId = ret.workflowId
        if( !workflowId )
            throw new AbortOperationException("Invalid Seqera Platform API response - Missing workflow Id")
        if( ret.message )
            log.warn(ret.message.toString())

        // Prepare to collect report paths if tower configuration has a 'reports' section
        reports.flowCreate(workflowId)
    }

    protected void setAuthToken(SimpleHttpClient client, String token) {
        // check for plain jwt token
        if( token.count('.')==2 ) {
            client.setBearerToken(token)
            return
        }

        // try checking personal access token
        try {
            final plain = new String(token.decodeBase64())
            final p = plain.indexOf('.')
            if( p!=-1 && new JsonSlurper().parseText(  plain.substring(0, p) )  ) {
                // ok this is bearer token
                client.setBearerToken(token)
                return
            }
        }
        catch ( Exception e ) {
            log.trace "Enable to set bearer token ~ Reason: $e.message"
        }

        // fallback on simple token
        client.setBasicToken(TOKEN_PREFIX + token)
    }

    protected Map makeCreateReq(Session session) {
        def result = new HashMap(5)
        result.sessionId = session.uniqueId.toString()
        result.runName = session.runName
        result.projectName = session.workflowMetadata.projectName
        result.repository = session.workflowMetadata.repository
        result.workflowId = env.get('TOWER_WORKFLOW_ID')
        result.instant = Instant.now().toEpochMilli()
        this.towerLaunch = result.workflowId != null
        return result
    }

    @Override
    void onProcessCreate(TaskProcessor process){
        log.trace "Creating process ${process.name}"
        if( !processNames.add(process.name) )
            throw new IllegalStateException("Process name `${process.name}` already used")
    }

    @Override
    void onFlowBegin() {
        // configure error retry
        httpClient.maxRetries = maxRetries
        httpClient.backOffBase = backOffBase
        httpClient.backOffDelay = backOffDelay

        final req = makeBeginReq(session)
        final resp = sendHttpMessage(urlTraceBegin, req, 'PUT')
        if( resp.error ) {
            log.debug """\
                Unexpected HTTP response
                - endpoint    : $urlTraceBegin
                - status code : $resp.code
                - response msg: $resp.cause
                """.stripIndent(true)
            throw new AbortOperationException(resp.message)
        }

        final payload = parseTowerResponse(resp)
        this.watchUrl = payload.watchUrl
        this.sender = Threads.start('Tower-thread', this.&sendTasks0)
        final msg = "Monitor the execution with Seqera Platform using this URL: ${watchUrl}"
        log.info(LoggerHelper.STICKY, msg)
    }

    String getAccessToken() {
        // when 'TOWER_WORKFLOW_ID' is provided in the env, it's a tower made launch
        // therefore the access token should only be taken from the env
        // otherwise check into the config file and fallback in the env
        def token = env.get('TOWER_WORKFLOW_ID')
                ? env.get('TOWER_ACCESS_TOKEN')
                : session.config.navigate('tower.accessToken', env.get('TOWER_ACCESS_TOKEN'))
        if( !token )
            throw new AbortOperationException("Missing personal access token -- Make sure there's a variable TOWER_ACCESS_TOKEN in your environment")
        return token
    }

    /**
     * Send an HTTP message when the workflow is completed.
     */
    @Override
    void onFlowComplete() {
        // submit the record
        events << new ProcessEvent(completed: true)
        // publish runtime reports
        reports.publishRuntimeReports()
        // wait the submission of pending events
        sender.join()
        // wait and flush reports content
        reports.flowComplete()
        // notify the workflow completion
        terminated = true
        final req = makeCompleteReq(session)
        final resp = sendHttpMessage(urlTraceComplete, req, 'PUT')
        logHttpResponse(urlTraceComplete, resp)
    }

    @Override
    void onProcessPending(TaskHandler handler, TraceRecord trace) {
        events << new ProcessEvent(trace: trace)
    }

    /**
     * Send an HTTP message when a process has been submitted
     *
     * @param handler A {@link TaskHandler} object representing the task submitted
     * @param trace A {@link TraceRecord} object holding the task metadata and runtime info
     */
    @Override
    void onProcessSubmit(TaskHandler handler, TraceRecord trace) {
        events << new ProcessEvent(trace: trace)
    }

    /**
     * Send an HTTP message, when a process has started
     *
     * @param handler A {@link TaskHandler} object representing the task started
     * @param trace A {@link TraceRecord} object holding the task metadata and runtime info
     */
    @Override
    void onProcessStart(TaskHandler handler, TraceRecord trace) {
        events << new ProcessEvent(trace: trace)
    }

    /**
     * Send an HTTP message, when a process completed
     *
     * @param handler A {@link TaskHandler} object representing the task completed
     * @param trace A {@link TraceRecord} object holding the task metadata and runtime info
     */
    @Override
    void onProcessComplete(TaskHandler handler, TraceRecord trace) {
        events << new ProcessEvent(trace: trace)

        synchronized (this) {
            aggregator.aggregate(trace)
        }

    }

    @Override
    void onProcessCached(TaskHandler handler, TraceRecord trace) {
        // event was triggered by a stored task, ignore it
        if( trace == null )
            return

        // add the cached task event
        events << new ProcessEvent(trace: trace)

        // remove the record from the current records
        synchronized (this) {
            aggregator.aggregate(trace)
        }
    }

    /**
     * Send an HTTP message, when a workflow has failed
     *
     * @param handler A {@link TaskHandler} object representing the task that caused the workflow execution to fail (it may be null)
     * @param trace A {@link TraceRecord} object holding the task metadata and runtime info (it may be null)
     */
    @Override
    void onFlowError(TaskHandler handler, TraceRecord trace) {
        events << new ProcessEvent(trace: trace)
    }

    /**
     * Update reports file when a file is published
     *
     * @param destination File path at `publishDir` of the published file.
     */
    @Override
    void onFilePublish(Path destination) {
        reports.filePublish(destination)
    }

    protected void refreshToken(String refresh) {
        log.debug "Token refresh request >> $refresh"
        final url = "$endpoint/oauth/access_token"
        httpClient.sendHttpMessage(
                url,
                method: 'POST',
                contentType: "application/x-www-form-urlencoded",
                body: "grant_type=refresh_token&refresh_token=${URLEncoder.encode(refresh, 'UTF-8')}" )

        final authCookie = httpClient.getCookie('JWT')
        final refreshCookie = httpClient.getCookie('JWT_REFRESH_TOKEN')

        // set the new bearer token
        if( authCookie?.value ) {
            log.trace "Updating http client bearer token=$authCookie.value"
            httpClient.setBearerToken(authCookie.value)
        }
        else {
            log.warn "Missing JWT cookie from refresh token response ~ $authCookie"
        }

        // set the new refresh token
        if( refreshCookie?.value ) {
            log.trace "Updating http client refresh token=$refreshCookie.value"
            refreshToken = refreshCookie.value
        }
        else {
            log.warn "Missing JWT_REFRESH_TOKEN cookie from refresh token response ~ $refreshCookie"
        }
    }

    /**
     * Little helper method that sends a HTTP POST message as JSON with
     * the current run status, ISO 8601 UTC timestamp, run name and the TraceRecord
     * object, if present.
     * @param event The current run status. One of {'started', 'process_submit', 'process_start',
     * 'process_complete', 'error', 'completed'}
     * @param payload An additional object to send. Must be of type TraceRecord or Manifest
     */
    protected Response sendHttpMessage(String url, Map payload, String method='POST'){

        int refreshTries=0
        final currentRefresh = refreshToken ?: env.get('TOWER_REFRESH_TOKEN')

        while ( true ) {
            // The actual HTTP request
            final String json = payload != null ? generator.toJson(payload) : null
            final String debug = json != null ? JsonOutput.prettyPrint(json).indent() : '-'
            log.trace "HTTP url=$url; payload:\n${debug}\n"
            try {
                if( refreshTries==1 ) {
                    refreshToken(currentRefresh)
                }

                httpClient.sendHttpMessage(url, json, method)
                return new Response(httpClient.responseCode, httpClient.getResponse())
            }
            catch( ConnectException e ) {
                String msg = "Unable to connect to Seqera Platform API: ${getHostUrl(url)}"
                return new Response(0, msg)
            }
            catch (IOException e) {
                int code = httpClient.responseCode
                if( code == 401 && ++refreshTries==1 && currentRefresh ) {
                    // when 401 Unauthorized error is returned - only the very first time -
                    // and a refresh token is available, make another iteration trying
                    // having refreshed the authorization token (see 'refreshToken' invocation above)
                    log.trace "Got 401 Unauthorized response ~ tries refreshing auth token"
                    continue
                }
                else {
                    log.trace("Got HTTP code $code - refreshTries=$refreshTries - currentRefresh=$currentRefresh", e)
                }

                String msg
                if( code == 401 ) {
                    msg = 'Unauthorized Seqera Platform API access -- Make sure you have specified the correct access token'
                }
                else {
                    msg = parseCause(httpClient.response) ?: "Unexpected response for request $url"
                }
                return new Response(code, msg, httpClient.response)
            }
        }
    }

    protected boolean isCliLogsEnabled() {
        return env.get('TOWER_ALLOW_NEXTFLOW_LOGS') == 'true'
    }

    protected String getOperationId() {
        if( !isCliLogsEnabled() )
            return null
        try {
            if( env.get('AWS_BATCH_JOB_ID') )
                return  "aws-batch::${env.get('AWS_BATCH_JOB_ID')}"
            else
                return "local-platform::${ProcessHelper.selfPid()}"
        }
        catch (Exception e) {
            log.warn "Unable to retrieve native environment operation id", e
            return null
        }
    }

    protected String getLogFile() {
        return isCliLogsEnabled() ? env.get('NXF_LOG_FILE') : null
    }

    protected String getOutFile() {
        return isCliLogsEnabled() ? env.get('NXF_OUT_FILE') : null
    }

    protected Map makeBeginReq(Session session) {
        def workflow = session.getWorkflowMetadata().toMap()
        workflow.params = session.getParams()
        workflow.id = getWorkflowId()
        workflow.remove('stats')

        // render as a string
        workflow.container = mapToString(workflow.container)
        workflow.configText = session.resolvedConfig
        // extra metadata
        workflow.operationId = getOperationId()
        workflow.logFile = getLogFile()
        workflow.outFile = getOutFile()

        def result = new LinkedHashMap(5)
        result.workflow = workflow
        result.processNames = new ArrayList(processNames)
        result.towerLaunch = towerLaunch
        result.instant = Instant.now().toEpochMilli()
        return result
    }

    protected Map makeCompleteReq(Session session) {
        def workflow = session.getWorkflowMetadata().toMap()
        workflow.params = session.getParams()
        workflow.id = getWorkflowId()
        // render as a string
        workflow.container = mapToString(workflow.container)
        workflow.configText = session.resolvedConfig
        // extra metadata
        workflow.operationId = getOperationId()
        workflow.logFile = getLogFile()
        workflow.outFile = getOutFile()

        def result = new LinkedHashMap(5)
        result.workflow = workflow
        result.metrics = getMetricsList()
        result.progress = getWorkflowProgress(false)
        result.instant = Instant.now().toEpochMilli()
        return result
    }

    protected Map makeHeartbeatReq() {
        def result = new HashMap(1)
        result.progress = getWorkflowProgress(true)
        result.instant = Instant.now().toEpochMilli()
        return result
    }

    protected String mapToString(def obj) {
        if( obj == null )
            return null
        if( obj instanceof CharSequence )
            return obj.toString()
        if( obj instanceof Map ) {
            // turn this off for multiple containers because the string representation is broken
            return null
        }
        throw new IllegalArgumentException("Illegal container attribute type: ${obj.getClass().getName()} = ${obj}" )
    }

    protected Map makeTaskMap0(TraceRecord trace) {
        Map<String,?> record = new LinkedHashMap<>(trace.store.size())
        for( Map.Entry<String,Object> entry : trace.store.entrySet() ) {
            def name = entry.key
            // remove '%' char from field prefix
            if( name.startsWith('%') )
                name = 'p' + name.substring(1)
            // normalise to camelCase
            name = underscoreToCamelCase(name)
            // put the value
            record.put(name, fixTaskField(name,entry.value))
        }

        // prevent invalid tag data
        if( record.tag!=null && !(record.tag instanceof CharSequence)) {
            final msg = "Invalid tag value for process: ${record.process} -- A string is expected instead of type: ${record.tag.getClass().getName()}; offending value=${record.tag}"
            log.warn1(msg, cacheKey: record.process)
            record.tag = null
        }

        // add transient fields
        record.executor = trace.getExecutorName()
        record.cloudZone = trace.getMachineInfo()?.zone
        record.machineType = trace.getMachineInfo()?.type
        record.priceModel = trace.getMachineInfo()?.priceModel?.toString()

        return record
    }


    static protected Object fixTaskField(String name, value) {
        if( TraceRecord.FIELDS[name] == 'date' )
            return value ? OffsetDateTime.ofInstant(Instant.ofEpochMilli(value as long), ZoneId.systemDefault()) : null
        else
            return value
    }

    protected Map makeTasksReq(Collection<TraceRecord> tasks) {

        def payload = new ArrayList(tasks.size())
        for( TraceRecord rec : tasks ) {
            payload << makeTaskMap0(rec)
        }

        final result = new LinkedHashMap(5)
        result.put('tasks', payload)
        result.put('progress', getWorkflowProgress(true))
        result.instant = Instant.now().toEpochMilli()
        return result
    }

    protected List getMetricsList() {
        return aggregator.computeSummaryList()
    }

    protected WorkflowProgress getWorkflowProgress(boolean quick) {
        def stats = quick ? session.getStatsObserver().getQuickStats() : session.getStatsObserver().getStats()
        new WorkflowProgress(stats)
    }

    /**
     * Little helper function that can be called for logging upon an incoming HTTP response
     */
    protected void logHttpResponse(String url, Response resp){
        if (resp.code >= 200 && resp.code < 300) {
            log.trace "Successfully send message to ${url} -- received status code ${resp.code}"
        }
        else {
            def cause = parseCause(resp.cause)
            def msg = """\
                Unexpected HTTP response.
                Failed to send message to ${endpoint} -- received 
                - status code : $resp.code
                - response msg: $resp.message
                """.stripIndent(true)
            // append separately otherwise formatting get broken
            msg += "- error cause : ${cause ?: '-'}"
            log.warn(msg)
        }
    }

    protected Map parseTowerResponse(Response resp) {
        if (resp.code >= 200 && resp.code < 300) {
            return (Map)new JsonSlurper().parseText(resp.message)
        }

        def cause = parseCause(resp.cause)

        def msg = """\
                Unexpected Seqera Platform API response
                - endpoint url: $endpoint
                - status code : $resp.code
                - response msg: ${resp.message} 
                """.stripIndent(true)
        // append separately otherwise formatting get broken
        msg += "- error cause : ${cause ?: '-'}"
        throw new Exception(msg)
    }

    protected String parseCause(String cause) {
        if( !cause )
            return null
        try {
            def map = (Map)new JsonSlurper().parseText(cause)
            return map.message
        }
        catch ( Exception ) {
            log.debug "Unable to parse error cause as JSON object: $cause"
            return cause
        }
    }

    protected String underscoreToCamelCase(String str) {
        if( !str.contains('_') )
            return str

        final words = str.tokenize('_')
        def result = words[0]
        for( int i=1; i<words.size(); i++ )
            result+=words[i].capitalize()

        return result
    }


    protected Map<String,Integer> loadSchema() {
        final props = new Properties()
        props.load(this.getClass().getResourceAsStream('/tower-schema.properties'))
        final result = new HashMap<String,Integer>(props.size())
        for( String key : props.keySet() ) {
            final value = props.getProperty(key)
            result.put( key, value ? value as Integer : null )
        }
        return result
    }

    protected void sendTasks0(dummy) {
        final tasks = new HashMap<TaskId, TraceRecord>(TASKS_PER_REQUEST)
        boolean complete = false
        long previous = System.currentTimeMillis()
        final long period = requestInterval.millis
        final long delay = period / 10 as long

        while( !complete ) {
            final ProcessEvent ev = events.poll(delay, TimeUnit.MILLISECONDS)
            // reconcile task events ie. send out only the last event
            if( ev ) {
                log.trace "Tower event=$ev"
                if( ev.trace )
                    tasks[ev.trace.taskId] = ev.trace
                if( ev.completed )
                    complete = true
            }

            // check if there's something to send
            final now = System.currentTimeMillis()
            final delta = now -previous

            if( !tasks ) {
                if( delta > aliveInterval.millis ) {
                    final req = makeHeartbeatReq()
                    final resp = sendHttpMessage(urlTraceHeartbeat, req, 'PUT')
                    logHttpResponse(urlTraceHeartbeat, resp)
                    previous = now
                }
                continue
            }

            if( delta > period || tasks.size() >= TASKS_PER_REQUEST || complete ) {
                // send
                final req = makeTasksReq(tasks.values())
                final resp = sendHttpMessage(urlTraceProgress, req, 'PUT')
                logHttpResponse(urlTraceProgress, resp)

                // clean up for next iteration
                previous = now
                tasks.clear()
            }
        }
    }

}
