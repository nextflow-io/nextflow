/*
 * Copyright 2022, Google Inc.
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

package nextflow.cloud.google.batch.client

import com.google.auth.oauth2.GoogleCredentials
import groovy.json.JsonOutput
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.cloud.google.batch.json.JsonHelper
import nextflow.cloud.google.batch.model.BatchJob
/**
 * Implements Google Batch HTTP client
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class BatchClient {

    protected String projectId
    protected String location
    protected GoogleCredentials credentials
    protected String BATCH_ENDPOINT = 'https://batch.googleapis.com/v1alpha1'
    protected String LOGS_ENDPOINT = 'https://logging.googleapis.com'

    BatchClient(BatchConfig config) {
        this.projectId = config.googleOpts.projectId
        this.location = config.googleOpts.location
        this.credentials = config.credentials
    }

    /** Only for testing - do not use */
    protected BatchClient() {}

    protected String getToken() {
        // only for debugging purpose
        if( System.getenv('_GOOGLE_ACCESS_TOKEN') ) {
            return System.getenv('_GOOGLE_ACCESS_TOKEN')
        }
        // use native API
        credentials.refreshIfExpired()
        final result = credentials.getAccessToken()?.getTokenValue()
        return result
    }

    BatchJsonResponse submitJob(String jobId, BatchJob request) {
        final action = "$BATCH_ENDPOINT/projects/${projectId}/locations/${location}/jobs?job_id=${jobId}"
        final resp = post(action, request.toJson())
        trace('POST', action, resp.text)
        return new BatchJsonResponse(resp.text)
    }

    BatchJsonResponse describeJob(String jobId) {
        final action = "$BATCH_ENDPOINT/projects/${projectId}/locations/$location/jobs/$jobId"
        final resp = get(action)
        trace('GET', action, resp.text)
        return new BatchJsonResponse(resp.text)
    }

    BatchJsonResponse listTasks(String jobId) {
        final action = "$BATCH_ENDPOINT/projects/${projectId}/locations/${location}/jobs/${jobId}/taskGroups/group0/tasks"
        final resp = get(action)
        trace('GET', action, resp.text)
        return new BatchJsonResponse(resp.text)
    }

    BatchJsonResponse deleteJobs(String jobId) {
        final action = "$BATCH_ENDPOINT/projects/${projectId}/locations/${location}/jobs/${jobId}"
        final resp = delete(action)
        trace('DELETE', action, resp.text)
        return new BatchJsonResponse(resp.text)
    }

    BatchJsonResponse logsEntries(String pageToken=null) {
        final action = "$LOGS_ENDPOINT/v2/entries:list"
        final req = [resourceNames: ["projects/$projectId"], filter:"logs/batch_task_logs", pageToken: pageToken]
        final resp = post(action, JsonHelper.toJson(req))
        trace('POST', action, resp.text)
        return new BatchJsonResponse(resp.text)
    }

    Map getJobStatus(String jobId) {
        final resp = describeJob(jobId)
        return (Map) resp.status
    }

    String getJobState(String jobId) {
        final status = getJobStatus(jobId)
        return status ? status.state : null
    }

    protected BatchApiResponse get(String path) {
        makeRequest('GET',path)
    }

    protected BatchApiResponse post(String path, String spec) {
        makeRequest('POST', path, spec)
    }

    protected BatchApiResponse delete(String path, String body=null) {
        makeRequest('DELETE', path, body)
    }

    /**
     * Makes a HTTP(S) request to the Batch service
     *
     * @param method The HTTP verb to use eg. {@code GET}, {@code POST}, etc
     * @param path The API action path
     * @param body The request payload
     * @return
     *      A two elements list in which the first entry is an integer representing the HTTP response code,
     *      the second element is the text (json) response
     */
    protected BatchApiResponse makeRequest(String method, String path, String body=null) throws BatchResponseException {

        final int maxAttempts = 5
        int attempt = 0

        while ( attempt < maxAttempts ) {
            attempt++

            try {
                return makeRequestCall( method, path, body )
            }
            catch ( SocketException e ) {
                log.warn "[GOOGLE BATCH] API request threw socket exception: $e.message for $method $path ${body ? '\n'+prettyPrint(body).indent() : ''}"
                if ( attempt < maxAttempts ) log.info( "[GOOGLE BATCH] Try API request again, remaining attempts: ${ maxAttempts - attempt }" )
                else throw e
                final long delay = (Math.pow(3, attempt - 1) as long) * 250
                sleep( delay )
            }
        }
    }


    private BatchApiResponse makeRequestCall(String method, String path, String body=null) throws BatchResponseException {
        final conn = createConnection0(path)
        conn.setRequestProperty("Content-Type", "application/json")
        conn.setRequestProperty("Authorization", "Bearer ${getToken()}")
        conn.setRequestProperty("X-Goog-User-Project", projectId)
        
        if( !method ) method = body ? 'POST' : 'GET'
        conn.setRequestMethod(method)
        log.trace "[GOOGLE BATCH] API request $method $path ${body ? '\n'+prettyPrint(body).indent() : ''}"

        if( body ) {
            conn.setDoOutput(true);
            conn.setDoInput(true);
            conn.getOutputStream() << body
            conn.getOutputStream().flush()
        }

        final code = conn.getResponseCode()
        final isError = code >= 400
        final stream = isError ? conn.getErrorStream() : conn.getInputStream()
        if( isError )
            throw new BatchResponseException("Request $method $path returned an error code=$code", stream)
        return new BatchApiResponse(code, stream)
    }

    protected HttpURLConnection createConnection0(String url) {
        new URL(url).openConnection() as HttpURLConnection
    }

    static private void trace(String method, String path, String text) {
        log.trace "[GOOGLE BATCH] API response $method $path \n${prettyPrint(text).indent()}"
    }

    static protected String prettyPrint(String json) {
        try {
            JsonOutput.prettyPrint(json)
        }
        catch( Exception e ) {
            return json
        }
    }

}
