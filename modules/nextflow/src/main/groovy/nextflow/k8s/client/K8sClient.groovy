/*
 * Copyright 2020-2022, Seqera Labs
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

package nextflow.k8s.client

import nextflow.exception.K8sOutOfCpuException
import nextflow.exception.K8sOutOfMemoryException

import javax.net.ssl.HostnameVerifier
import javax.net.ssl.HttpsURLConnection
import javax.net.ssl.SSLContext
import javax.net.ssl.SSLSession
import javax.net.ssl.TrustManager
import javax.net.ssl.TrustManagerFactory
import javax.net.ssl.X509TrustManager
import java.nio.file.Path
import java.security.KeyStore
import java.security.SecureRandom
import java.security.cert.CertificateFactory
import java.security.cert.X509Certificate

import groovy.json.JsonOutput
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.exception.NodeTerminationException
import nextflow.exception.ProcessFailedException
import org.yaml.snakeyaml.Yaml
/**
 * Kubernetes API client
 *
 * Tip: use the following command to find out your kubernetes master node
 *   kubectl cluster-info
 *
 * See
 *   https://v1-8.docs.kubernetes.io/docs/api-reference/v1.8/#-strong-api-overview-strong-
 *
 * Useful cheatsheet
 *   https://kubernetes.io/docs/reference/kubectl/cheatsheet/
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class K8sClient {

    protected ClientConfig config

    private TrustManager[] trustManagers

    private HostnameVerifier hostnameVerifier

    K8sClient() {
        this(new ClientConfig())
    }

    ClientConfig getConfig() { config }

    /**
     * Creates a kubernetes client using the configuration setting provided by the specified
     * {@link ConfigDiscovery} instance
     *
     * @param config
     */
    K8sClient(ClientConfig config) {
        this.config = config
        setupSslCert()
    }


    protected setupSslCert() {

        if( !config.verifySsl ) {
            // -- no SSL is required - use fake trust manager
            final trustAll = new X509TrustManager() {
                @Override X509Certificate[] getAcceptedIssuers() { return null }
                @Override void checkClientTrusted(X509Certificate[] certs, String authType) { }
                @Override void checkServerTrusted(X509Certificate[] certs, String authType) { }
            }

            trustManagers = [trustAll] as TrustManager[]
            hostnameVerifier = new HostnameVerifier() {
                @Override boolean verify(String hostname, SSLSession session) { return true }
            }

        }
        else if ( config.sslCert != null) {
            char[] password = null
            final factory = CertificateFactory.getInstance("X.509");
            final authority = new ByteArrayInputStream(config.sslCert)
            final certificates = factory.generateCertificates(authority)
            if (certificates.isEmpty()) {
                throw new IllegalArgumentException("Trusted certificates set cannot be empty");
            }

            final keyStore = KeyStore.getInstance(KeyStore.getDefaultType());
            keyStore.load(null, password);
            certificates.eachWithIndex{ cert, index ->
                String alias = "ca$index"
                keyStore.setCertificateEntry(alias, cert);
            }

            final trustManagerFactory = TrustManagerFactory.getInstance(TrustManagerFactory.getDefaultAlgorithm());
            trustManagerFactory.init(keyStore);
            trustManagers = trustManagerFactory.getTrustManagers();
        }
    }

    K8sResponseJson secretesList() {
        final action = "/api/v1/namespaces/$config.namespace/secrets"
        final resp = get(action)
        trace('GET', action, resp.text)
        new K8sResponseJson(resp.text)
    }

    K8sResponseJson secretDescribe(String name) {
        assert name
        final action = "/api/v1/namespaces/$config.namespace/secrets/$name"
        final resp = get(action)
        trace('GET', action, resp.text)
        new K8sResponseJson(resp.text)
    }

    /**
     * Create a pod
     *
     * See
     *  https://v1-8.docs.kubernetes.io/docs/api-reference/v1.8/#create-55
     *  https://v1-8.docs.kubernetes.io/docs/api-reference/v1.8/#pod-v1-core
     *
     * @param spec
     * @return
     */
    K8sResponseJson podCreate(String req) {
        assert req
        final action = "/api/v1/namespaces/$config.namespace/pods"
        final resp = post(action, req)
        trace('POST', action, resp.text)
        return new K8sResponseJson(resp.text)
    }

    K8sResponseJson podCreate(Map req, Path saveYamlPath=null) {

        if( saveYamlPath ) try {
            saveYamlPath.text = new Yaml().dump(req).toString()
        }
        catch( Exception e ) {
            log.debug "WARN: unable to save request yaml -- cause: ${e.message ?: e}"
        }

        podCreate(JsonOutput.toJson(req))
    }

    /**
     * Delete a pod
     *
     * See
     *   https://v1-8.docs.kubernetes.io/docs/api-reference/v1.8/#delete-58
     *
     * @param name
     * @return
     */
    K8sResponseJson podDelete(String name) {
        assert name
        final action = "/api/v1/namespaces/$config.namespace/pods/$name"
        final resp = delete(action)
        trace('DELETE', action, resp.text)
        new K8sResponseJson(resp.text)
    }

    /**
     * Create a job
     *
     * See
     *  https://v1-8.docs.kubernetes.io/docs/api-reference/v1.8/#create-55
     *  https://v1-8.docs.kubernetes.io/docs/api-reference/v1.8/#job-v1-batch
     *
     * @param spec
     * @return
     */
    K8sResponseJson jobCreate(String req) {
        assert req
        final action = "/apis/batch/v1/namespaces/$config.namespace/jobs"
        final resp = post(action, req)
        trace('POST', action, resp.text)
        return new K8sResponseJson(resp.text)
    }

    K8sResponseJson jobCreate(Map req, Path saveYamlPath=null) {

        if( saveYamlPath ) try {
            saveYamlPath.text = new Yaml().dump(req).toString()
        }
        catch( Exception e ) {
            log.debug "WARN: unable to save request yaml -- cause: ${e.message ?: e}"
        }

        jobCreate(JsonOutput.toJson(req))
    }

    /**
     * Delete a job
     *
     * See
     *   https://v1-8.docs.kubernetes.io/docs/api-reference/v1.8/#delete-58
     *
     * @param name
     * @return
     */

    K8sResponseJson jobDelete(String name) {
        assert name

        // get podList of a job
        final action = "/api/v1/namespaces/$config.namespace/pods?labelSelector=job-name=$name"
        final resp = get(action)
        trace('GET', action, resp.text)
        final podList = new K8sResponseJson(resp.text)

        // delete all pods in a job
        if (podList.kind == "PodList") { 
            for (item in podList.items) {
                try {
                   podDelete(((item as Map).metadata as Map).name as String)
                }
                catch(K8sResponseException err) {
                   if( err.response.code == 404 )
                       log.debug("Unable to delete Pod for job $name, pod already gone")
                   else
                       throw err
                }
            }
        }

        // delete job
        final action1 = "/apis/batch/v1/namespaces/$config.namespace/jobs/$name"
        final resp1 = delete(action1)
        trace('DELETE', action1, resp1.text)
        new K8sResponseJson(resp1.text)
    }

    /*
     * https://v1-8.docs.kubernetes.io/docs/api-reference/v1.8/#list-62
     * https://v1-8.docs.kubernetes.io/docs/api-reference/v1.8/#list-all-namespaces-63
     */

    K8sResponseJson jobList(boolean allNamespaces=false) {
        final String action = allNamespaces ? "jobs" : "namespaces/$config.namespace/jobs"
        final resp = get("/apis/batch/v1/$action")
        trace('GET', action, resp.text)
        new K8sResponseJson(resp.text)
    }

    /*
     * https://v1-8.docs.kubernetes.io/docs/api-reference/v1.8/#list-62
     * https://v1-8.docs.kubernetes.io/docs/api-reference/v1.8/#list-all-namespaces-63
     */
    K8sResponseJson podList(boolean allNamespaces=false) {
        final String action = allNamespaces ? "pods" : "namespaces/$config.namespace/pods"
        final resp = get("/api/v1/$action")
        trace('GET', action, resp.text)
        new K8sResponseJson(resp.text)
    }

    /*
     * https://v1-8.docs.kubernetes.io/docs/api-reference/v1.8/#read-status-69
     */

    // converts the name of job to the name of the latest created pod
    String findPodNameForJob(String name){
        final action = "/api/v1/namespaces/$config.namespace/pods?labelSelector=job-name=$name"
        final resp = get(action)
        trace('GET', action, resp.text)
        final podList = new K8sResponseJson(resp.text)

        String podName

        // find latest created pod
        if (podList.kind == "PodList") {
            final pods = podList.items
            String latestPod = "0000-00-00T00:00:00Z"

            for (item in pods) {
                final podMetadata = (Map) (item as Map).metadata
                final podTimestamp = podMetadata.creationTimestamp

                if (podTimestamp.toString() > latestPod) {
                    latestPod = podTimestamp
                    podName = podMetadata.name
                }
            }
        }
        return podName
    }

    K8sResponseJson jobStatus(String name) {
        assert name
        final action = "/apis/batch/v1/namespaces/$config.namespace/jobs/$name/status"
        final resp = get(action)
        trace('GET', action, resp.text)
        return new K8sResponseJson(resp.text)
    }

    K8sResponseJson podStatus(String name) {
        assert name
        final action = "/api/v1/namespaces/$config.namespace/pods/$name/status"
        final resp = get(action)
        trace('GET', action, resp.text)
        return new K8sResponseJson(resp.text)
    }

    /*
     * https://v1-8.docs.kubernetes.io/docs/api-reference/v1.8/#read-status-69
     */

    protected K8sResponseJson podStatus0(String name) {
        try {
            return podStatus(name)
        }
        catch (K8sResponseException err) {
            if( err.response.code == 404 && isKindPods(err.response)  ) {
                // this may happen when K8s node is shutdown and the pod is evicted
                // therefore process exception is thrown so that the failure
                // can be managed by the nextflow as re-triable execution
                throw new NodeTerminationException("Unable to find pod $name - The pod may be evicted by a node shutdown event")
            }
            throw err
        }
    }

    protected boolean isKindPods(K8sResponseJson resp) {
        if( resp.details instanceof Map ) {
            final details = (Map) resp.details
            return details.kind == 'pods'
        }
        return false
    }

    String getNodeOfPod(String podName){
        assert podName
        final K8sResponseJson resp = podStatus0(podName)
        (resp?.spec as Map)?.nodeName as String
    }

    /**
     * Get pod current state object
     *
     * @param podName The pod name
     * @return
     *      A {@link Map} representing the container state object as shown below
     *      <code>
     *       {
     *                "terminated": {
     *                    "exitCode": 127,
     *                    "reason": "ContainerCannotRun",
     *                    "message": "OCI runtime create failed: container_linux.go:296: starting container process caused \"exec: \\\"bash\\\": executable file not found in $PATH\": unknown",
     *                    "startedAt": "2018-01-12T22:04:25Z",
     *                    "finishedAt": "2018-01-12T22:04:25Z",
     *                    "containerID": "docker://730ef2e05be72ffc354f2682b4e8300610812137b9037b726c21e5c4e41b6dda"
     *                }
     *      </code>
     *      See the following link for details https://kubernetes.io/docs/reference/generated/kubernetes-api/v1.10/#containerstate-v1-core
     *      An empty map is return if the pod is a `Pending` status and the container state is not
     *      yet available
     *
     *
     */
    Map jobState( String jobName ) {
        assert jobName
        final podName = findPodNameForJob(jobName)
        if( podName ) {
            try {
                return podState(podName)
            } 
            /* pod might be deleted by control plane just after findPodNameForJob() call
             * so try fallback to jobState
             */   
            catch (NodeTerminationException err) {
                log.warn1("Job $jobName's Pod not found, probably cleaned by controlplane")
                return jobStateFallback0(jobName)           
            }
        }
        else {
            return jobStateFallback0(jobName)
        }
    }

    protected Map jobStateFallback0(String jobName) {
        final K8sResponseJson jobResp = jobStatus(jobName)
        final jobStatus = jobResp.status as Map
        if( jobStatus?.succeeded == 1 && jobStatus.conditions instanceof List ) {
            final allConditions = jobStatus.conditions as List<Map>
            final cond = allConditions.find { cond -> cond.type == 'Complete' }

            if( cond?.status == 'True' ) {
                log.warn1("Job $jobName already completed and Pod is gone")
                final dummyPodStatus = [
                        terminated: [
                                exitcode: 0,
                                reason: "Completed",
                                startedAt: jobStatus.startTime,
                                finishedAt: jobStatus.completionTime,
                        ]
                ]
                return dummyPodStatus
            } else {
                throw new ProcessFailedException("K8s Job $jobName succeeded but does not have Complete status. $allConditions")
            }
        }

        if( jobStatus?.failed && (int)(jobStatus.failed) > 0 ) {
            String message = 'unknown'
            if( jobStatus.conditions instanceof List ) {
                final allConditions = jobStatus.conditions as List<Map>
                final cond = allConditions.find { cond -> cond.type == 'Failed' }
                message = cond?.message
            }
            throw new ProcessFailedException("K8s Job $jobName execution failed: $message")
        }

        log.warn1("K8s Job $jobName does not have pod - Not yet scheduled?")
        return Collections.emptyMap()
    }

    /**
     * Get pod current state object
     *
     * @param podName The pod name
     * @return
     *      A {@link Map} representing the container state object as shown below
     *      <code>
     *       {
     *                "terminated": {
     *                    "exitCode": 127,
     *                    "reason": "ContainerCannotRun",
     *                    "message": "OCI runtime create failed: container_linux.go:296: starting container process caused \"exec: \\\"bash\\\": executable file not found in $PATH\": unknown",
     *                    "startedAt": "2018-01-12T22:04:25Z",
     *                    "finishedAt": "2018-01-12T22:04:25Z",
     *                    "containerID": "docker://730ef2e05be72ffc354f2682b4e8300610812137b9037b726c21e5c4e41b6dda"
     *                }
     *      </code>
     *      See the following link for details https://kubernetes.io/docs/reference/generated/kubernetes-api/v1.10/#containerstate-v1-core
     *      An empty map is return if the pod is a `Pending` status and the container state is not
     *      yet available
     *
     *
     */
    Map podState( String podName ) {
        assert podName

        final K8sResponseJson resp = podStatus0(podName)
        final status = resp.status as Map
        final containerStatuses = status?.containerStatuses as List<Map>

        if( containerStatuses?.size()>0 ) {
            final container = containerStatuses.get(0)
            // note: when the pod is created by a Job submission
            // the `podName` does not match the container name because it
            // contains a suffix random generated by K8s pod scheduler
            if( !container.name || !podName.startsWith(container.name.toString()) )
                throw new K8sResponseException("K8s invalid status for pod: ${podName} (unexpected container name: ${container.name})", resp)

            if( !container.state )
                throw new K8sResponseException("K8s invalid status for pod: ${podName} (missing state object)", resp)

            final state = container.state as Map
            if( state.waiting instanceof Map ) {
                def waiting = state.waiting as Map
                checkInvalidWaitingState(waiting, resp)
            }
            return state
        }

        if( status?.phase == 'Pending' ){
            if( status.conditions instanceof List ) {
                final allConditions = status.conditions as List<Map>
                final cond = allConditions.find { cond -> cond.type == 'PodScheduled' }
                if( cond.reason == 'Unschedulable' ) {
                    def message = "K8s pod cannot be scheduled"
                    if( cond.message ) message += " -- $cond.message"
                    //def cause = new K8sResponseException(resp)
                    log.warn1(message)
                }
            }
            // undetermined status -- return an empty response
            return Collections.emptyMap()
        }

        if( status?.phase == 'Failed' ) {
            def msg = "K8s pod '$podName' execution failed"
            if( status.reason ) msg += " - reason: ${status.reason}"
            if( status.message ) msg += " - message: ${status.message}"
            switch ( status.reason ) {
                case 'OutOfcpu':    throw new K8sOutOfCpuException(msg)
                case 'OutOfmemory': throw new K8sOutOfMemoryException(msg)
                case 'Shutdown':    throw new NodeTerminationException(msg)
                default:            throw new ProcessFailedException(msg)
            }
        }

        throw new K8sResponseException("K8s undetermined status conditions for pod $podName", resp)
    }

    protected void checkInvalidWaitingState( Map waiting, K8sResponseJson resp ) {
        if( waiting.reason == 'ErrImagePull' || waiting.reason == 'ImagePullBackOff') {
            def message = "K8s pod image cannot be pulled"
            if( waiting.message ) message += " -- $waiting.message"
            final cause = new K8sResponseException(resp)
            throw new PodUnschedulableException(message, cause)
        }
        if( waiting.reason == 'CreateContainerConfigError' ) {
            def message = "K8s pod configuration failed"
            if( waiting.message ) message += " -- $waiting.message"
            final cause = new K8sResponseException(resp)
            throw new PodUnschedulableException(message, cause)
        }
        if( waiting.reason =~ /.+Error$/ ) {
            def message = "K8s pod waiting on unknown error state"
            if( waiting.message ) message += " -- $waiting.message"
            final cause = new K8sResponseException(resp)
            throw new PodUnschedulableException(message, cause)
        }
        final status = resp.status as Map
        if( status?.phase == 'Failed' ) {
            def message = "K8s pod in Failed state"
            final cause = new K8sResponseException(resp)
            throw new PodUnschedulableException(message, cause)
        }
    }

    /*
     * https://v1-8.docs.kubernetes.io/docs/api-reference/v1.8/#read-log
     */
    InputStream jobLog(String name) {
        jobLog( Collections.emptyMap(), name )
    }

    InputStream jobLog(Map params, String name) {
        assert name
        // -- compose the request action uri
        String podName = findPodNameForJob(name)
        def action = "/api/v1/namespaces/$config.namespace/pods/$podName/log"
        int count=0
        for( String key : (params.keySet()) ) {
            action += "${count++==0 ? '?' : '&'}${key}=${params.get(key)}"
        }
        // -- submit request
        def resp = get(action)
        resp.stream
    }

    /*
     * https://v1-8.docs.kubernetes.io/docs/api-reference/v1.8/#read-log
     */
    InputStream podLog(String name) {
        podLog( Collections.emptyMap(), name )
    }

    InputStream podLog(Map params, String name) {
        assert name
        // -- compose the request action uri
        def action = "/api/v1/namespaces/$config.namespace/pods/$name/log"
        int count=0
        for( String key : (params.keySet()) ) {
            action += "${count++==0 ? '?' : '&'}${key}=${params.get(key)}"
        }
        // -- submit request
        def resp = get(action)
        resp.stream
    }

    protected K8sResponseApi post(String path, String spec) {
        makeRequest('POST', path, spec)
    }

    protected K8sResponseApi delete(String path, String body=null) {
        makeRequest('DELETE', path, body)
    }

    protected HttpURLConnection createConnection0(String url) {
        new URL(url).openConnection() as HttpURLConnection
    }

    protected void setupHttpsConn( HttpsURLConnection conn ) {
        if (config.httpReadTimeout != null) {
            conn.setReadTimeout(config.httpReadTimeout.toMillis() as int)
        }
        if (config.httpConnectTimeout != null) {
            conn.setConnectTimeout(config.httpConnectTimeout.toMillis() as int)
        }
        if (config.keyManagers != null || trustManagers != null) {
            SSLContext sslContext = SSLContext.getInstance("TLS");
            sslContext.init(config.keyManagers, trustManagers, new SecureRandom());
            conn.setSSLSocketFactory(sslContext.getSocketFactory())
        }

        if( hostnameVerifier )
            conn.setHostnameVerifier(hostnameVerifier)
    }

    /**
     * Makes a HTTP(S) request the kubernetes master
     *
     * @param method The HTTP verb to use eg. {@code GET}, {@code POST}, etc
     * @param path The API action path
     * @param body The request payload
     * @return
     *      A two elements list in which the first entry is an integer representing the HTTP response code,
     *      the second element is the text (json) response
     */
    protected K8sResponseApi makeRequest(String method, String path, String body=null) throws K8sResponseException {

        final int maxRetries = config.maxErrorRetry
        int attempt = 0

        while ( true ) {
            try {
                return makeRequestCall( method, path, body )
            } catch ( K8sResponseException | SocketException | SocketTimeoutException e ) {
                if ( e instanceof K8sResponseException && e.response.code != 500 )
                    throw e
                if ( ++attempt > maxRetries )
                    throw e
                log.debug "[K8s] API request threw socket exception: $e.message for $method $path - Retrying request (attempt=$attempt)"
                final long delay = (Math.pow(3, attempt - 1) as long) * 250
                sleep( delay )
            }
        }
    }


    private K8sResponseApi makeRequestCall(String method, String path, String body=null) throws K8sResponseException {
        assert config.server, 'Missing Kubernetes server name'
        assert path.startsWith('/'), 'Kubernetes API request path must starts with a `/` character'

        final prefix = config.server.contains("://") ? config.server : "https://$config.server"
        final conn = createConnection0(prefix + path)
        conn.setRequestProperty("Content-Type", "application/json")
        if( config.token ) {
            conn.setRequestProperty("Authorization", "Bearer $config.token")
        }

        if( conn instanceof HttpsURLConnection ) {
            setupHttpsConn(conn)
        }

        if( !method ) method = body ? 'POST' : 'GET'
        conn.setRequestMethod(method)
        log.trace "[K8s] API request $method $path ${body ? '\n'+prettyPrint(body).indent() : ''}"

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
            throw new K8sResponseException("Request $method $path returned an error code=$code", stream)
        return new K8sResponseApi(code, stream)
    }

    static private void trace(String method, String path, String text) {
        log.trace "[K8s] API response $method $path \n${prettyPrint(text).indent()}"
    }

    protected K8sResponseApi get(String path) {
        makeRequest('GET',path)
    }

    static protected String prettyPrint(String json) {
        try {
            JsonOutput.prettyPrint(json)
        }
        catch( Exception e ) {
            return json
        }
    }

    K8sResponseJson configCreate(String name, Map data) {

        final spec = [
                apiVersion: 'v1',
                kind: 'ConfigMap',
                metadata: [ name: name, namespace: config.namespace ],
                data: data
        ]

        configCreate0(spec)
    }

    protected K8sResponseJson configCreate0(Map spec) {
        final action = "/api/v1/namespaces/${config.namespace}/configmaps"
        final body = JsonOutput.toJson(spec)
        def resp = post(action, body)
        trace('POST', action, resp.text)
        return new K8sResponseJson(resp.text)
    }


    K8sResponseJson configDelete(String name) {
        final action = "/api/v1/namespaces/${config.namespace}/configmaps/$name"
        def resp = delete(action)
        trace('DELETE', action, resp.text)
        return new K8sResponseJson(resp.text)
    }

    K8sResponseJson configDeleteAll() {
        final action = "/api/v1/namespaces/${config.namespace}/configmaps"
        def resp = delete(action)
        trace('DELETE', action, resp.text)
        return new K8sResponseJson(resp.text)
    }


    K8sResponseJson volumeClaimRead(String name) {
        final action = "/api/v1/namespaces/${config.namespace}/persistentvolumeclaims/${name}"
        def resp = get(action)
        trace('GET', action, resp.text)
        return new K8sResponseJson(resp.text)
    }


}

