/*
 * Copyright (c) 2019, Seqera Labs.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * This Source Code Form is "Incompatible With Secondary Licenses", as
 * defined by the Mozilla Public License, v. 2.0.
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
import groovy.transform.PackageScope
import groovy.transform.ToString
import groovy.transform.TupleConstructor
import groovy.util.logging.Slf4j
import nextflow.Const
import nextflow.NextflowMeta
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
import nextflow.util.SimpleHttpClient
/**
 * Send out messages via HTTP to a configured URL on different workflow
 * execution events.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class TowerObserver implements TraceObserver {

    static final String DEF_ENDPOINT_URL = 'https://api.tower.nf'

    static private final int TASKS_PER_REQUEST = 100

    static private final Duration REQUEST_INTERVAL = Duration.of('1 sec')

    static private final Duration ALIVE_INTERVAL = Duration.of('1 min')

    static private final String TOKEN_PREFIX = '@token:'

    @TupleConstructor
    static class Response {
        final int code
        final String message
        boolean isError() { code < 200 || code >= 300 }
    }

    @ToString(includeNames = true)
    class ProcessEvent {
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
    protected SimpleHttpClient httpClient

    private JsonGenerator generator = createJsonGeneratorForPayloads()

    private String workflowId

    private String watchUrl

    private String endpoint

    private String urlTraceWorkflow

    private String urlTraceTask

    private String urlTraceAlive

    private String urlTraceInit

    private ResourcesAggregator aggregator

    private Map<String,String> env = System.getenv()

    private LinkedBlockingQueue<ProcessEvent> events = new LinkedBlockingQueue()

    private Thread sender

    private Duration requestInterval = REQUEST_INTERVAL

    private Duration aliveInterval = ALIVE_INTERVAL

    private LinkedHashSet<String> processNames = new LinkedHashSet<>(20)

    private boolean terminated

    /**
     * Constructor that consumes a URL and creates
     * a basic HTTP client.
     * @param endpoint The target address for sending messages to
     */
    TowerObserver(String endpoint) {
        this.endpoint = checkUrl(endpoint)
        this.urlTraceTask = this.endpoint + '/trace/task'
        this.urlTraceWorkflow = this.endpoint + '/trace/workflow'
        this.urlTraceAlive = this.endpoint + '/trace/alive'
        this.urlTraceInit = this.endpoint + '/trace/init'
    }

    /**
     * only for testing purpose -- do not use
     */
    protected TowerObserver() {

    }

    boolean enableMetrics() { true }

    String getEndpoint() { endpoint }

    void setAliveInterval(Duration d) {
        this.aliveInterval = d
    }

    void setRequestInterval(Duration d ) {
        this.requestInterval = d
    }

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
        if( url =~ "^(https|http)://[-a-zA-Z0-9+&@#/%?=~_|!:,.;]*[-a-zA-Z0-9+&@#/%=~_|]" ) {
            while( url.endsWith('/') )
                url = url[0..-2]
            return url
        }
        throw new IllegalArgumentException("Only http or https are supported protocols -- The given URL was: ${url}")
    }

    protected String getHostUrl(String endpoint) {
        def url = new URL(endpoint)
        return "${url.protocol}://${url.authority}"
    }

    /**
     * On workflow start, submit a message with some basic
     * information, like Id, activity and an ISO 8601 formatted
     * timestamp.
     * @param session The current Nextflow session object
     */
    @Override
    void onFlowInit(Session session) {
        log.debug "Creating Tower observer -- endpoint=$endpoint; requestInterval=$requestInterval; aliveInterval=$aliveInterval"
        
        this.session = session
        this.aggregator = new ResourcesAggregator(session)
        this.runName = session.getRunName()
        this.runId = session.getUniqueId()
        this.httpClient = new SimpleHttpClient().setAuthToken(TOKEN_PREFIX + getAccessToken())

        // send hello to verify auth
        final msg = [sessionId: session.uniqueId]
        final resp = sendHttpMessage(urlTraceInit, msg)
        if( resp.error )
            throw new AbortOperationException(resp.message)
        final ret = parseTowerResponse(resp)
        this.workflowId = ret.workflowId
        if( !workflowId )
            throw new AbortOperationException("Invalid Tower response")
    }

    @Override
    void onProcessCreate(TaskProcessor process){
        log.debug "Creating process ${process.name}"
        if( !processNames.add(process.name) )
            throw new IllegalStateException("Process name `${process.name}` already used")
    }

    @Override
    void onFlowBegin() {
        final req = makeWorkflowReq(session)
        final resp = sendHttpMessage(urlTraceWorkflow, req)
        if( resp.error )
            throw new AbortOperationException(resp.message)

        final payload = parseTowerResponse(resp)
        this.watchUrl = payload.watchUrl
        this.sender = Thread.start('Tower-thread', this.&sendTasks0)
        final msg = "Monitor the execution with Nextflow Tower using this url ${watchUrl}"
        log.info(LoggerHelper.STICKY, msg)
    }

    protected String getAccessToken() {
        // access token
        def token = session.config.navigate('tower.accessToken')
        if( !token )
            token = env.get('TOWER_ACCESS_TOKEN')
        if( !token )
            throw new AbortOperationException("Missing Nextflow Tower access token -- Make sure there's a variable TOWER_ACCESS_TOKEN in your environment")
        return token
    }

    /**
     * Send an HTTP message when the workflow is completed.
     */
    @Override
    void onFlowComplete() {
        // submit the record
        events << new ProcessEvent(completed: true)
        // wait the submission of pending events
        sender.join()
        // notify the workflow completion
        terminated = true
        def resp = sendHttpMessage(urlTraceWorkflow, makeWorkflowReq(session))
        logHttpResponse(urlTraceWorkflow, resp)
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

    protected Response sendHttpGet(String url) {
        sendHttpMessage(url, null, 'GET')
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
        // The actual HTTP request
        final String json = payload != null ? generator.toJson(payload) : null
        final String debug = json != null ? JsonOutput.prettyPrint(json).indent() : '-'
        log.debug "HTTP url=$url; payload:\n${debug}\n"
        try {
            httpClient.sendHttpMessage(url, json, method)
            return new Response(httpClient.responseCode, httpClient.getResponse())
        }
        catch( ConnectException e ) {
            String msg = "Unable to connect Tower host: ${getHostUrl(url)}"
            return new Response(0, msg)
        }
        catch (IOException e) {
            int code = httpClient.responseCode
            String msg = ( code == 401
                            ? 'Unauthorized Tower access -- Make sure you have specified the correct access token'
                            : "Unexpected response code $code for request $url"  )
            return new Response(code, msg)
        }
    }

    protected Map makeWorkflowReq(Session session) {
        def workflow = session.getWorkflowMetadata().toMap()
        workflow.params = session.getParams()
        workflow.id = workflowId
        if( !terminated )
            workflow.remove('stats')
        // render as a string
        workflow.container = mapToString(workflow.container)
        workflow.configText = session.resolvedConfig

        def result = new LinkedHashMap(5)
        result.workflow = workflow
        if( terminated ) {
            result.metrics = getMetricsList()
        }
        else {
            result.processNames = new ArrayList(processNames)
        }
        return result
    }

    protected String mapToString(def obj) {
        if( obj == null )
            return null
        if( obj instanceof CharSequence )
            return obj.toString()
        if( obj instanceof Map ) {
            def map = obj as Map
            return map.collect { k,v -> "$k:$v" }.join(',')
        }
        throw new IllegalArgumentException("Illegal container attribut type: ${obj.getClass().getName()} = ${obj}" )
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
        result.put('workflowId', workflowId)
        return result
    }

    protected List getMetricsList() {
        return aggregator.computeSummaryList()
    }

    /**
     * Little helper function that can be called for logging upon an incoming HTTP response
     */
    protected void logHttpResponse(String url, Response resp){
        if (resp.code >= 200 && resp.code < 300) {
            log.trace "Successfully send message to ${url} -- received status code ${resp.code}"
        }
        else {
            def msg = """\
                Unexpected HTTP response.
                Failed to send message to ${endpoint} -- received 
                - status code : $resp.code    
                - response msg: $resp.message   
                """.stripIndent()
            log.warn(msg)
        }
    }

    protected Map parseTowerResponse(Response resp) {
        if (resp.code >= 200 && resp.code < 300) {
            return (Map)new JsonSlurper().parseText(resp.message)
        }

        def msg = """\
                Unexpected Tower response
                - endpoint url: $endpoint
                - status code : $resp.code    
                - response msg: ${resp.message}  
                """.stripIndent()

        throw new Exception(msg)
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

    @PackageScope
    static JsonGenerator createJsonGeneratorForPayloads() {
        new JsonGenerator.Options()
                .addConverter(Path) { Path p, String key -> p.toUriString() }
                .addConverter(Duration) { Duration d, String key -> d.durationInMillis }
                .addConverter(NextflowMeta) { NextflowMeta m, String key -> m.toJsonMap() }
                .addConverter(OffsetDateTime) { it.toString() }
                .dateFormat(Const.ISO_8601_DATETIME_FORMAT).timezone("UTC")
                .build()
    }


    protected void sendTasks0(dummy) {
        final tasks = new HashMap<TaskId, TraceRecord>(TASKS_PER_REQUEST)
        boolean complete = false
        long previous = System.currentTimeMillis()
        final long period = requestInterval.millis
        final long delay = period / 10 as long

        while( !complete ) {
            final ev = events.poll(delay, TimeUnit.MILLISECONDS)
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
                    final req = new HashMap(1); req.workflowId = workflowId
                    final resp = sendHttpMessage(urlTraceAlive, req)
                    logHttpResponse(urlTraceAlive, resp)
                    previous = now
                }
                continue
            }

            if( delta > period || tasks.size() >= TASKS_PER_REQUEST || complete ) {
                // send
                final req = makeTasksReq(tasks.values())
                final resp = sendHttpMessage(urlTraceTask, req)
                logHttpResponse(urlTraceTask, resp)

                // clean up for next iteration
                previous = now
                tasks.clear()
            }
        }
    }

}
