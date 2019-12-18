/*
 * Copyright 2018, University of Tübingen, Quantitative Biology Center (QBiC)
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

package nextflow.trace

import groovy.json.JsonGenerator
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import groovyx.gpars.agent.Agent
import nextflow.Const
import nextflow.NextflowMeta
import nextflow.Session
import nextflow.processor.TaskHandler
import nextflow.script.ScriptBinding.ParamsMap
import nextflow.script.WorkflowMetadata
import nextflow.util.Duration
import nextflow.util.SimpleHttpClient

import java.nio.file.Path

/**
 * Send out messages via HTTP to a configured URL on different workflow
 * execution events.
 *
 * @author Sven Fillinger <sven.fillinger@qbic.uni-tuebingen.de>
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class WebLogObserver implements TraceObserver{

    private Session session

    private static final TimeZone UTC = TimeZone.getTimeZone("UTC")

    /**
     * Workflow identifier, will be taken from the Session() object later
     */
    private String runName

    /**
     * Store the sessions unique ID for downstream reference purposes
     */
    private String runId

    /**
     * The default url is localhost
     */
    static public String DEF_URL = 'http://localhost'

    /**
     * Simple http client object that will send out messages
     */
    private SimpleHttpClient httpClient = new SimpleHttpClient()

    /**
     * An agent for the http request in an own thread
     */
    private Agent<WebLogObserver> webLogAgent

    /**
     * Json generator for weblog payloads
     */
    private JsonGenerator generator

    private String endpoint

    /**
     * Constructor that consumes a URL and creates
     * a basic HTTP client.
     * @param url The target address for sending messages to
     */
    WebLogObserver(String url) {
        this.endpoint = checkUrl(url)
        this.webLogAgent = new Agent<>(this)
        this.generator = createJsonGeneratorForPayloads()
    }

    /**
     * only for testing purpose -- do not use
     */
    protected WebLogObserver() {

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
            return url
        }
        throw new IllegalArgumentException("Only http or https are supported protocols -- The given URL was: ${url}")
    }

    /**
     * On workflow start, submit a message with some basic
     * information, like Id, activity and an ISO 8601 formatted
     * timestamp.
     * @param session The current Nextflow session object
     */
    @Override
    void onFlowCreate(Session session) {
        this.session = session
        runName = session.getRunName()
        runId = session.getUniqueId()

        asyncHttpMessage("started", createFlowPayloadFromSession(session))
    }

    /**
     * Send an HTTP message when the workflow is completed.
     */
    @Override
    void onFlowComplete() {
        asyncHttpMessage("completed", createFlowPayloadFromSession(this.session))
        webLogAgent.await()
    }

    /**
     * Send an HTTP message when a process has been submitted
     *
     * @param handler A {@link TaskHandler} object representing the task submitted
     * @param trace A {@link TraceRecord} object holding the task metadata and runtime info
     */
    @Override
    void onProcessSubmit(TaskHandler handler, TraceRecord trace) {
        asyncHttpMessage("process_submitted", trace)
    }

    /**
     * Send an HTTP message, when a process has started
     *
     * @param handler A {@link TaskHandler} object representing the task started
     * @param trace A {@link TraceRecord} object holding the task metadata and runtime info
     */
    @Override
    void onProcessStart(TaskHandler handler, TraceRecord trace) {
        asyncHttpMessage("process_started", trace)
    }

    /**
     * Send an HTTP message, when a process completed
     *
     * @param handler A {@link TaskHandler} object representing the task completed
     * @param trace A {@link TraceRecord} object holding the task metadata and runtime info
     */
    @Override
    void onProcessComplete(TaskHandler handler, TraceRecord trace) {
        asyncHttpMessage("process_completed", trace)
    }

    /**
     * Send an HTTP message, when a workflow has failed
     *
     * @param handler A {@link TaskHandler} object representing the task that caused the workflow execution to fail (it may be null)
     * @param trace A {@link TraceRecord} object holding the task metadata and runtime info (it may be null)
     */
    @Override
    void onFlowError(TaskHandler handler, TraceRecord trace) {
        asyncHttpMessage("error", trace)
    }

    /**
     * Little helper method that sends a HTTP POST message as JSON with
     * the current run status, ISO 8601 UTC timestamp, run name and the TraceRecord
     * object, if present.
     * @param event The current run status. One of {'started', 'process_submit', 'process_start',
     * 'process_complete', 'error', 'completed'}
     * @param payload An additional object to send. Must be of type TraceRecord or Manifest
     */
    protected void sendHttpMessage(String event, Object payload = null){

        // Set the message info
        final time = new Date().format(Const.ISO_8601_DATETIME_FORMAT, UTC)

        final message = new HashMap(4)
        message.runName = runName
        message.runId = runId
        message.event = event
        message.utcTime = time

        if (payload instanceof TraceRecord)
                message.trace = (payload as TraceRecord).store
        else if (payload instanceof FlowPayload)
                message.metadata = payload
        else if (payload != null)
            throw new IllegalArgumentException("Only TraceRecord and Manifest class types are supported: [${payload.getClass().getName()}] $payload")

        // The actual HTTP request
        httpClient.sendHttpMessage(endpoint, generator.toJson(message))
        logHttpResponse()
    }

    protected static FlowPayload createFlowPayloadFromSession(Session session) {
        def params = session.binding.getProperty('params') as ParamsMap
        def workflow = session.getWorkflowMetadata()
        new FlowPayload(params, workflow)
    }

    /**
     * Asynchronous HTTP POST request wrapper.
     * @param event The workflow run status
     * @param payload An additional object to send. Must be of type TraceRecord or Manifest
     */
    protected void asyncHttpMessage(String event, Object payload = null){
        webLogAgent.send{sendHttpMessage(event, payload)}
    }

    /**
     * Little helper function that can be called for logging upon an incoming HTTP response
     */
    protected void logHttpResponse(){
        def statusCode = httpClient.getResponseCode()
        if (statusCode >= 200 && statusCode < 300) {
            log.debug "Successfully send message to ${endpoint} -- received status code ${statusCode}."
        } else {
            def msg = """\
                Unexpected HTTP response.
                Failed to send message to ${endpoint} -- received 
                - status code : $statusCode    
                - response msg: ${httpClient.getResponse()}  
                """.stripIndent()
            log.debug msg
        }
    }

    private static JsonGenerator createJsonGeneratorForPayloads() {
        new JsonGenerator.Options()
                .addConverter(Path) { Path p, String key -> p.toUriString() }
                .addConverter(Duration) { Duration d, String key -> d.durationInMillis }
                .addConverter(NextflowMeta) { meta, key -> meta.toJsonMap() }
                .dateFormat(Const.ISO_8601_DATETIME_FORMAT).timezone("UTC")
                .build()
    }

    private static class FlowPayload {

        final ParamsMap parameters

        final WorkflowMetadata workflow

        FlowPayload(ParamsMap params, WorkflowMetadata workflow ) {
            this.parameters = params
            this.workflow = workflow
        }
    }
}
