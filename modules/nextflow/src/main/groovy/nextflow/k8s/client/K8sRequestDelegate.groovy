/*
 * Copyright 2013-2023, Seqera Labs
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

import javax.net.ssl.HostnameVerifier
import javax.net.ssl.HttpsURLConnection
import javax.net.ssl.SSLContext
import javax.net.ssl.SSLSession
import javax.net.ssl.TrustManager
import javax.net.ssl.TrustManagerFactory
import javax.net.ssl.X509TrustManager
import java.security.KeyStore
import java.security.SecureRandom
import java.security.cert.CertificateFactory
import java.security.cert.X509Certificate

import groovy.json.JsonOutput
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
/**
 * Implements helper methods for classes that use the Kubernetes REST API.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class K8sRequestDelegate {

    private ClientConfig config

    private TrustManager[] trustManagers

    private HostnameVerifier hostnameVerifier

    ClientConfig getConfig() { config }

    K8sRequestDelegate() {
        this(new ClientConfig())
    }

    K8sRequestDelegate(ClientConfig config) {
        this.config = config
        setupSslCert()
    }

    protected void setupSslCert() {

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

    /**
     * Make a HTTP(S) request to the Kubernetes cluster.
     *
     * @param method The HTTP verb to use eg. {@code GET}, {@code POST}, etc
     * @param path The API action path
     * @param body The request payload
     * @return
     *      A two elements list in which the first entry is an integer representing the HTTP response code,
     *      the second element is the text (json) response
     */
    K8sResponseApi makeRequest(String method, String path, String body=null) throws K8sResponseException {

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

    void trace(String method, String path, String text) {
        log.trace "[K8s] API response $method $path \n${prettyPrint(text).indent()}"
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

