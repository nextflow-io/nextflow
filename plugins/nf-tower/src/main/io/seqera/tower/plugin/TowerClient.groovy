/*
 * Copyright 2013-2026, Seqera Labs
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

package io.seqera.tower.plugin

import java.net.http.HttpClient
import java.net.http.HttpRequest
import java.time.Instant
import java.time.OffsetDateTime
import java.time.ZoneId
import java.util.concurrent.ConcurrentHashMap
import java.util.concurrent.LinkedBlockingQueue
import java.util.concurrent.TimeUnit

import groovy.json.JsonGenerator
import groovy.json.JsonOutput
import groovy.json.JsonSlurper
import groovy.transform.CompileStatic
import groovy.transform.ToString
import groovy.transform.TupleConstructor
import groovy.util.logging.Slf4j
import io.seqera.http.HxClient
import io.seqera.util.trace.TraceUtils
import nextflow.BuildInfo
import nextflow.Session
import nextflow.container.resolver.ContainerMeta
import nextflow.exception.AbortOperationException
import nextflow.processor.TaskHandler
import nextflow.processor.TaskId
import nextflow.processor.TaskProcessor
import nextflow.trace.ResourcesAggregator
import nextflow.trace.TraceObserverV2
import nextflow.trace.TraceRecord
import nextflow.trace.event.FilePublishEvent
import nextflow.trace.event.TaskEvent
import nextflow.util.Duration
import nextflow.util.LoggerHelper
import nextflow.util.ProcessHelper
import nextflow.util.TestOnly
import nextflow.util.Threads
/**
 * Perform HTTP call to Seqera platform.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class TowerClient {

    static final public String DEF_ENDPOINT_URL = 'https://api.cloud.seqera.io'

    static private final String TOKEN_PREFIX = '@token:'

    @TupleConstructor
    static class Response {
        final int code
        final String message
        final String cause
        boolean isError() { code < 200 || code >= 300 }
    }

    private HxClient httpClient

    private JsonGenerator generator

    private String endpoint

    protected Map<String,String> env = System.getenv()

    private Map<String,Integer> schema = Collections.emptyMap()

    private String accessToken

    private TowerRetryPolicy retryPolicy
    private Duration readTimeout = TowerConfig.DEFAULT_READ_TIMEOUT

    private Duration connectTimeout = TowerConfig.DEFAULT_CONNECT_TIMEOUT

    TowerClient(TowerConfig config, Map env) {
        this.endpoint = checkUrl(config.endpoint)
        this.accessToken = config.accessToken
        this.retryPolicy = config.retryPolicy
        this.readTimeout = config.httpReadTimeout
        this.connectTimeout = config.httpConnectTimeout
        this.schema = loadSchema()
        this.generator = TowerJsonGenerator.create(schema)
        if (env)
            this.env = env
        initHttpClient()
    }

    TowerClient withRequestTimeout(Duration duration) {
        this.readTimeout = duration
        return this
    }

    @TestOnly
    protected TowerClient() {
        this.generator = TowerJsonGenerator.create(Collections.EMPTY_MAP)
    }

    String getEndpoint() { endpoint }

    /**
     * Check the URL and create an HttpPost() object. If a invalid i.e. protocol is used,
     * the constructor will raise an exception.
     *
     * The RegEx was taken and adapted from http://urlregex.com
     *
     * @param url String with target URL
     * @return The requested url or the default url, if invalid
     */
    protected String checkUrl(String url) {
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

    Map traceCreate(Map req, String workspaceId){
        return sendAndProcessRequest( getUrlTraceCreate(workspaceId), req, 'POST')
    }

    Map traceBegin(Map req, String workspaceId, String workflowId){
        return sendAndProcessRequest( getUrlTraceBegin(workspaceId, workflowId), req, 'POST')
    }

    void traceComplete(Map req, String workspaceId, String workflowId) {
        final url = getUrlTraceComplete(workspaceId, workflowId)
        final resp = sendHttpMessage(url, req, 'PUT')
        logHttpResponse(url, resp)
    }

    void traceHeartbeat(Map req, String workspaceId, String workflowId) {
        final url = getUrlTraceHeartbeat( workspaceId, workflowId)
        final resp = sendHttpMessage(url, req, 'PUT')
        logHttpResponse(url, resp)
    }

    void traceProgress(Map req, String workspaceId, String workflowId) {
        final url = getUrlTraceProgress( workspaceId, workflowId )
        final resp = sendHttpMessage(url, req, 'PUT')
        logHttpResponse(url, resp)
    }

    protected Map sendAndProcessRequest(String url, Map req, String method){
        final resp = sendHttpMessage(url, req, method)
        if( resp.error ) {
            log.debug """\
                Unexpected HTTP response
                - endpoint    : $url
                - status code : $resp.code
                - response msg: $resp.cause
                """.stripIndent(true)
            throw new AbortOperationException(resp.message)
        }
        return parseTowerResponse(resp)
    }

    protected String getUrlTraceCreate(String workspaceId) {
        def result = this.endpoint + '/trace/create'
        if( workspaceId )
            result += "?workspaceId=$workspaceId"
        return result
    }

    protected String getUrlTraceBegin(String workspaceId, String workflowId) {
        def result = "$endpoint/trace/$workflowId/begin"
        if( workspaceId )
            result += "?workspaceId=$workspaceId"
        return result
    }

    protected String getUrlTraceComplete(String workspaceId, String workflowId) {
        def result = "$endpoint/trace/$workflowId/complete"
        if( workspaceId )
            result += "?workspaceId=$workspaceId"
        return result
    }

    protected String getUrlTraceHeartbeat(String workspaceId, String workflowId) {
        def result = "$endpoint/trace/$workflowId/heartbeat"
        if( workspaceId )
            result += "?workspaceId=$workspaceId"
        return result
    }

    protected String getUrlTraceProgress(String workspaceId, String workflowId) {
        def result = "$endpoint/trace/$workflowId/progress"
        if( workspaceId )
            result += "?workspaceId=$workspaceId"
        return result
    }

    protected void initHttpClient() {
        final builder = HxClient.newBuilder()
        // auth settings
        setupClientAuth(builder, getAccessToken())
        // retry settings
        this.httpClient = builder
            .retryConfig(this.retryPolicy)
            .followRedirects(HttpClient.Redirect.NORMAL)
            .version(HttpClient.Version.HTTP_1_1)
            .connectTimeout(java.time.Duration.ofMillis(connectTimeout.millis))
            .build()
    }

    protected void setupClientAuth(HxClient.Builder config, String token) {
        // check for plain jwt token
        final refreshToken = env.get('TOWER_REFRESH_TOKEN')
        final refreshUrl = refreshToken ? "$endpoint/oauth/access_token" : null
        if( token.count('.')==2 ) {
            config.bearerToken(token)
            config.refreshToken(refreshToken)
            config.refreshTokenUrl(refreshUrl)
            config.refreshCookiePolicy(CookiePolicy.ACCEPT_ALL)
            return
        }

        // try checking personal access token
        try {
            final plain = new String(token.decodeBase64())
            final p = plain.indexOf('.')
            if( p!=-1 && new JsonSlurper().parseText(  plain.substring(0, p) )  ) {
                // ok this is bearer token
                config.bearerToken(token)
                // setup the refresh
                config.refreshToken(refreshToken)
                config.refreshTokenUrl(refreshUrl)
                config.refreshCookiePolicy(CookiePolicy.ACCEPT_ALL)
                return
            }
        }
        catch ( Exception e ) {
            log.trace "Enable to set bearer token ~ Reason: $e.message"
        }

        // fallback on simple token
        config.basicAuth(TOKEN_PREFIX + token)
    }

    String getAccessToken() {
        if( !accessToken )
            throw new AbortOperationException("Missing Seqera Platform access token -- Make sure there's a variable TOWER_ACCESS_TOKEN in your environment")
        return accessToken
    }

    /**
     * Send an HTTP message, when a process has started
     *
     * @param handler A {@link TaskHandler} object representing the task started
     * @param trace A {@link TraceRecord} object holding the task metadata and runtime info
     */

    /**
     * Little helper method that sends a HTTP POST message as JSON with
     * the current run status, ISO 8601 UTC timestamp, run name and the TraceRecord
     * object, if present.
     * @param event The current run status. One of {'started', 'process_submit', 'process_start',
     * 'process_complete', 'error', 'completed'}
     * @param payload An additional object to send. Must be of type TraceRecord or Manifest
     */
    protected Response sendHttpMessage(String url, Map payload, String method='POST') {

        // The actual HTTP request
        final String json = payload != null ? generator.toJson(payload) : null
        final String debug = json != null ? JsonOutput.prettyPrint(json).indent() : '-'
        log.trace "HTTP url=$url; payload:\n${debug}\n"
        try {
            final resp = httpClient.sendAsString(makeRequest(url, json, method))
            final status = resp.statusCode()
            if( status == 401 ) {
                final msg = 'Unauthorized Seqera Platform API access -- Make sure you have specified the correct access token'
                return new Response(status, msg)
            }
            if( status>=400 ) {
                final msg = parseCause(resp?.body()) ?: "Unexpected response for request $url"
                return new Response(status, msg as String)
            }
            else {
                return new Response(status, resp.body())
            }
        }
        catch( IOException e ) {
            String msg = "Unable to connect to Seqera Platform API: ${getHostUrl(url)}"
            return new Response(0, msg)
        }
    }

    Response sendApiRequest(String url, Map payload=null, String method='GET') {
        sendHttpMessage(url, payload, method)
    }

    protected HttpRequest makeRequest(String url, String payload, String verb) {
        final builder = HttpRequest.newBuilder(URI.create(url))
            .header('User-Agent', "Nextflow/$BuildInfo.version")
            .header('Traceparent', TraceUtils.rndTrace())
            .timeout(java.time.Duration.ofMillis(readTimeout.millis))

        if( verb == 'GET' )
            return builder.GET().build()

        assert payload, "Tower request cannot be empty"
        builder.header('Content-Type', 'application/json; charset=utf-8')

        if( verb == 'PUT' )
            return builder.PUT(HttpRequest.BodyPublishers.ofString(payload)).build()

        if( verb == 'POST' )
            return builder.POST(HttpRequest.BodyPublishers.ofString(payload)).build()

        throw new IllegalArgumentException("Unsupported HTTP verb: $verb")
    }

    /**
     * Little helper function that can be called for logging upon an incoming HTTP response
     */
    protected void logHttpResponse(String url, Response resp) {
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

}
