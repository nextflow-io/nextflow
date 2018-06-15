/*
 * Copyright (c) 2018, University of Tübingen, Quantitative Biology Center (QBiC).
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
 */

package nextflow.trace

import groovy.json.JsonOutput
import groovy.json.JsonSlurper
import groovy.util.logging.Slf4j
import groovyx.gpars.agent.Agent

import nextflow.Session
import nextflow.script.WorkflowMetadata
import nextflow.util.SimpleHttpClient
import nextflow.processor.TaskHandler

/**
 * Send out messages via HTTP to a configured URL on different workflow
 * execution events.
 *
 * @author Sven Fillinger <sven.fillinger@qbic.uni-tuebingen.de>
 */
@Slf4j
class WebLogObserver implements TraceObserver{

    /**
     * Workflow identifier, will be taken from the Session() object later
     */
    private String runName

    /**
     * Store the sessions unique ID for downstream reference purposes
     */
    private String runId

    /**
     * Not a HTTP header, but a message header with some general workflow info
     */
    private static JsonSlurper JSLURPER = new JsonSlurper()

    /**
     * The default url is localhost
     */
    static String DEF_URL = 'http://localhost'

    /**
     * Contains server request response
     */
    String response = ""

    /**
     * Simple http client object that will send out messages
     */
    private SimpleHttpClient httpClient = new SimpleHttpClient()

    /**
     * Holds information of workflow metadata
     */
    private WorkflowMetadata workflowMetadata

    /**
     * An agent for the http request in an own thread
     */
    private Agent<WebLogObserver> webLogAgent

    /**
     * Constructor that consumes a URL and creates
     * a basic HTTP client.
     * @param url The target address for sending messages to
     */
    WebLogObserver(String url) {
        this.httpClient.setUrl(checkUrl(url))
        this.runName = ""
        this.runId = ""
        this.webLogAgent = new Agent<>(this)
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
    static String checkUrl(String url){
        try {
            assert url=~ "^(https|http)://[-a-zA-Z0-9+&@#/%?=~_|!:,.;]*[-a-zA-Z0-9+&@#/%=~_|]"
        } catch (AssertionError e){
            throw new IllegalArgumentException("Only http or https are supported protocols. The given URL was ${url}.")
        }
        return url
    }

    /**
     * On workflow start, submit a message with some basic
     * information, like Id, activity and an ISO 8601 formatted
     * timestamp.
     * @param session The current Nextflow session object
     */
    @Override
    void onFlowStart(Session session) {
        // This is either set by the user or via Nexflows name generator
        runName = session.getRunName()
        runId = session.getUniqueId()
        asyncHttpMessage("started")
    }

    /**
     * Send an HTTP message when the workflow is completed.
     */
    @Override
    void onFlowComplete() {
        asyncHttpMessage("completed")
    }

    /**
     * Send an HTTP message when a process has been submitted
     * @param handler
     * @param trace
     */
    @Override
    void onProcessSubmit(TaskHandler handler, TraceRecord trace) {
        asyncHttpMessage("process_submitted", trace)
    }

    /**
     * Send an HTTP message, when a process has started
     * @param handler
     * @param trace
     */
    @Override
    void onProcessStart(TaskHandler handler, TraceRecord trace) {
        asyncHttpMessage("process_started", trace)
    }

    /**
     * Send an HTTP message, when a process completed
     * @param handler
     * @param trace
     */
    @Override
    void onProcessComplete(TaskHandler handler, TraceRecord trace) {
        asyncHttpMessage("process_completed", trace)
    }

    /**
     * Send an HTTP message, when a workflow has failed
     * @param handler
     * @param trace
     */
    @Override
    void onFlowError(TaskHandler handler, TraceRecord trace) {
        asyncHttpMessage("error", trace)
    }

    /**
     * Little helper method that sends a HTTP POST message as JSON with
     * the current run status, ISO 8601 UTC timestamp, run name and the TraceRecord
     * object, if present.
     * @param runStatus The current run status. One of {'started', 'process_submit', 'process_start',
     * 'process_complete', 'error', 'completed'}
     * @param trace A TraceRecord object that contains current process information
     */
    protected void sendHttpMessage(String runStatus, TraceRecord trace = null){

        if (!this.httpClient.getUrl()){
            this.response = "No message send, url was not specified or formatted wrong."
            log.error(this.response)
            return
        }

        // Set the message info
        def time = new Date().format("yyyy-MM-dd'T'HH:mm:ss'Z'", TimeZone.getTimeZone("UTC"))
        def messageJson = JSLURPER.parseText(
                """{ "runName": "$runName", "runId": "$runId", "runStatus":"$runStatus", "utcTime":"$time" }""")

        // Append the trace object if present
        if (trace)
            messageJson["trace"] = JSLURPER.parseText(JsonOutput.toJson(trace.store))

        // The actual HTTP request
        httpClient.sendHttpMessage(JsonOutput.toJson(messageJson))
        logHttpResponse()

        this.response = httpClient.getResponse()
    }

    /**
     * Asynchronous HTTP POST request wrapper.
     * @param runStatus The workflow run status
     * @param trace A TraceRecord object with workflow information
     * @return A Java string, that contains the HTTP request response
     */
    protected void asyncHttpMessage(String runStatus, TraceRecord trace = null){
        webLogAgent.send{sendHttpMessage(runStatus, trace)}
    }

    /**
     * Little helper function that can be called for logging upon an incoming HTTP response
     */
    private void logHttpResponse(){
        def statusCode = httpClient.getResponseCode()
        if (statusCode == 200)
            log.debug "Successfully send message to ${httpClient.getUrl()}, received status code 200."
        else {
            log.debug "Failed to send message to ${httpClient.getUrl()}, received status code $statusCode"
            log.debug httpClient.getResponse()
        }
    }

}
