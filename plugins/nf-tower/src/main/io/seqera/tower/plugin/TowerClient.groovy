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

import groovy.json.JsonGenerator
import groovy.json.JsonOutput
import groovy.json.JsonSlurper
import groovy.transform.CompileStatic
import groovy.transform.TupleConstructor
import groovy.util.logging.Slf4j
import io.seqera.http.HxClient
import io.seqera.tower.plugin.exception.ForbiddenException
import io.seqera.tower.plugin.exception.NotFoundException
import io.seqera.tower.plugin.exception.UnauthorizedException
import io.seqera.util.trace.TraceUtils
import nextflow.BuildInfo
import nextflow.SysEnv
import nextflow.exception.AbortRunException
import nextflow.util.Duration
import nextflow.util.TestOnly
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

    private Map<String,Integer> schema = Collections.emptyMap()

    private String accessToken

    private TowerRetryPolicy retryPolicy

    private Duration readTimeout = TowerConfig.DEFAULT_READ_TIMEOUT

    private Duration connectTimeout = TowerConfig.DEFAULT_CONNECT_TIMEOUT

    TowerClient(TowerConfig config) {
        this.endpoint = checkUrl(config.endpoint)
        this.accessToken = config.accessToken
        this.retryPolicy = config.retryPolicy
        this.readTimeout = config.httpReadTimeout
        this.connectTimeout = config.httpConnectTimeout
        this.schema = loadSchema()
        this.generator = TowerJsonGenerator.create(schema)
        initHttpClient()
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
        return sendAndProcessRequest( getUrlTraceBegin(workspaceId, workflowId), req, 'PUT')
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
        if( resp.error ) {
            final message =  """\
                Unexpected HTTP response
                - endpoint    : $url
                - status code : $resp.code
                - response msg: $resp.message
                """.stripIndent(true)
            throw new AbortRunException(message)
        }
    }

    protected Map sendAndProcessRequest(String url, Map req, String method){
        final resp = sendHttpMessage(url, req, method)
        if( resp.error ) {
            final message =  """\
                Unexpected HTTP response
                - endpoint    : $url
                - status code : $resp.code
                - response msg: $resp.message
                """.stripIndent(true)
            throw new AbortRunException(message)
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
        final refreshToken = SysEnv.get('TOWER_REFRESH_TOKEN')
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
            throw new AbortRunException("Missing Seqera Platform access token -- Make sure there's a variable TOWER_ACCESS_TOKEN in your environment")
        return accessToken
    }

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

    /**
     * Send a GET request and return the response body as a streaming {@link InputStream}
     * instead of buffering the entire response into a {@link String}.
     * Uses {@code HxClient.sendAsStream()} which goes through the same retry and
     * auth chain as {@code sendAsString()}.
     *
     * Status codes are checked before returning — on error the stream is closed and
     * the same exceptions as {@link #checkResponse} are thrown.
     *
     * @param url the full API URL to GET
     * @return an InputStream over the response body
     * @throws UnauthorizedException on 401
     * @throws ForbiddenException on 403
     * @throws NotFoundException on 404
     */
    InputStream sendStreamingRequest(String url) {
        log.trace "HTTP streaming GET url=$url"
        final req = makeRequest(url, null, 'GET')
        final resp = httpClient.sendAsStream(req)
        final status = resp.statusCode()
        if( status >= 200 && status < 300 )
            return resp.body()
        // Error — close the stream and throw
        resp.body()?.close()
        if( status == 401 )
            throw new UnauthorizedException("Seqera authentication failed — check tower.accessToken or TOWER_ACCESS_TOKEN")
        if( status == 403 )
            throw new ForbiddenException("Forbidden — check permissions")
        if( status == 404 )
            throw new NotFoundException("Resource $url not found")
        throw new IOException("Seqera API error: HTTP ${status} for ${url}")
    }

    protected HttpRequest makeRequest(String url, String payload, String verb) {
        final builder = HttpRequest.newBuilder(URI.create(url))
            .header('User-Agent', "Nextflow/$BuildInfo.version")
            .header('Traceparent', TraceUtils.rndTrace())
            .timeout(java.time.Duration.ofMillis(readTimeout.millis))

        if( verb == 'GET' )
            return builder.GET().build()

        if( verb == 'DELETE' )
            return builder.DELETE().build()

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
        throw new AbortRunException(msg)
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

    String buildUrl( String path, Map queryParams) {
        def url = new StringBuilder(endpoint)
        if( !path.startsWith('/') ) {
            url.append('/')
        }
        url.append(path)

        if( queryParams && !queryParams.isEmpty() ) {
            url.append('?')
            url.append(queryParams.collect { k, v -> "${URLEncoder.encode(k.toString(), 'UTF-8')}=${URLEncoder.encode(v.toString(), 'UTF-8')}" }.join('&'))
        }

        return url.toString()
    }

    Map apiGet(String path, Map queryParams = [:]) {
        final url = buildUrl( path, queryParams)

        final response = sendApiRequest(url)
        checkResponse(response, url)
        return new JsonSlurper().parseText(response.message) as Map
    }

    Map apiPost(String path, Map queryParams, Map payload) {
        final url = buildUrl( path, queryParams)
        final response= sendApiRequest(url, payload, 'POST')
        checkResponse(response, url)
        return response.message ? new JsonSlurper().parseText(response.message) as Map : [:]
    }

    /**
     * @return current user info (id, userName, etc.) from GET /user-info
     */
    Map<String, Object> getUserInfo() {
        final json = apiGet("/user-info")
        return json.user as Map<String, Object>
    }

    /**
     * Calls the Seqera Platform to retrieve the workflow information.
     *
     * @param workflowId Id of the workflow
     * @return Map containing workflow information
     * @throws RuntimeException if the API call fails
     */
    Map getWorkflowDetails( String workflowId, Map queryParams = [:]) {
        final json = apiGet("/workflow/${workflowId}", queryParams)
        return json.workflow as Map
    }


    List<Map> listUserWorkspacesAndOrgs(String userId) {
        final json = apiGet("/user/${userId}/workspaces")
        return json.orgsAndWorkspaces as List<Map>
    }

    private static void checkResponse(Response resp, String url) {
        if (!resp.error) return
        final code = resp.code
        if (code == 401)
            throw new UnauthorizedException("Seqera authentication failed — check tower.accessToken or TOWER_ACCESS_TOKEN")
        if (code == 403)
            throw new ForbiddenException("Forbidden — check permissions")
        if (code == 404)
            throw new NotFoundException("Resource $url not found")
        throw new Exception("Seqera API error: HTTP ${code} for ${url}${resp.message ? ' - ' + resp.message :''}")
    }

    /**
     * Calls the Seqera Platform to retrieve the user's workspaces information
     * and select the one matching with the workspace Id.
     *
     * @param userId Id of the workspace user
     * @param workspaceId Id of the workspace
     * @return Map containing workspace information
     * @throws RuntimeException if the API call fails
     */
    Map getUserWorkspaceDetails( String userId, String workspaceId) {
        if( !userId || !workspaceId ) {
            return null
        }
        try {
            final orgsAndWorkspaces = listUserWorkspacesAndOrgs(userId)

            final workspace = orgsAndWorkspaces.find { ((Map) it).workspaceId?.toString() == workspaceId }
            if( workspace ) {
                final ws = workspace as Map
                return [
                    orgName          : ws.orgName,
                    workspaceId      : ws.workspaceId,
                    workspaceName    : ws.workspaceName,
                    workspaceFullName: ws.workspaceFullName,
                    roles            : ws.roles
                ]
            }

            return null
        } catch( Exception e ) {
            log.debug("Failed to get workspace details for workspace ${workspaceId}: ${e.message}", e)
            return null
        }
    }


}
