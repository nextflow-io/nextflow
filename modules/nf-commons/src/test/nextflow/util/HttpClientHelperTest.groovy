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

package nextflow.util

import java.net.http.HttpClient
import java.net.http.HttpRequest
import java.net.http.HttpResponse
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths
import java.security.KeyStore

import javax.net.ssl.KeyManagerFactory
import javax.net.ssl.SSLContext
import javax.net.ssl.TrustManagerFactory

import com.sun.net.httpserver.HttpServer
import com.sun.net.httpserver.HttpsConfigurator
import com.sun.net.httpserver.HttpsServer
import spock.lang.Shared
import spock.lang.Specification
import spock.util.environment.RestoreSystemProperties
import test.EnvHelper
import test.MockAuthProxyServer

/**
 * Verify HTTP clients configured via {@link HttpClientHelper#applyProxy} authenticate
 * to a forward proxy requiring Basic authentication, for both plain HTTP requests and
 * HTTPS CONNECT tunnelling.
 *
 * Note: the HTTPS tunnelling test requires the test JVM to be launched with
 * {@code -Djdk.http.auth.tunneling.disabledSchemes=} (see the test task in the
 * module build file) because by default the JDK strips Basic credentials from
 * the proxy CONNECT request. At runtime the property is set by the Nextflow
 * launcher when the proxy configuration provides credentials.
 *
 * @author Phil Ewels <phil.ewels@seqera.io>
 */
class HttpClientHelperTest extends Specification {

    static final Map<String,String> NO_PROXY_VARS = [
            HTTP_PROXY: null, http_proxy: null,
            HTTPS_PROXY: null, https_proxy: null,
            ALL_PROXY: null, all_proxy: null,
            NO_PROXY: null, no_proxy: null ]

    @Shared
    Path tempDir

    def setupSpec() {
        tempDir = Files.createTempDirectory('proxy-test')
    }

    def cleanupSpec() {
        tempDir?.deleteDir()
    }

    private static Map<String,String> proxyEnv(MockAuthProxyServer proxy, String user='foo', String password='secret') {
        final result = new HashMap<String,String>(NO_PROXY_VARS)
        result.HTTP_PROXY = "http://${user}:${password}@${proxy.host}:${proxy.port}".toString()
        result.HTTPS_PROXY = result.HTTP_PROXY
        result.NO_PROXY = ''
        return result
    }

    @RestoreSystemProperties
    def 'should authenticate a plain http request against the proxy'() {
        given:
        def proxy = new MockAuthProxyServer('foo', 'secret').start()
        proxy.responseBody = 'hello from proxy'

        when:
        HttpResponse<String> resp = null
        EnvHelper.withEnv(proxyEnv(proxy)) {
            final client = HttpClientHelper
                    .applyProxy(HttpClient.newBuilder().version(HttpClient.Version.HTTP_1_1))
                    .build()
            final request = HttpRequest.newBuilder()
                    .uri(new URI('http://server.proxy-test.internal/hello'))
                    .GET()
                    .build()
            resp = client.send(request, HttpResponse.BodyHandlers.ofString())
        }

        then: 'the request is answered by the proxy after the 407 challenge'
        resp.statusCode() == 200
        resp.body() == 'hello from proxy'
        and: 'the proxy challenged the client'
        proxy.unauthorizedCount.get() >= 1
        proxy.proxiedCount.get() >= 1

        cleanup:
        proxy?.close()
    }

    @RestoreSystemProperties
    def 'should fail when credentials are wrong'() {
        given:
        def proxy = new MockAuthProxyServer('foo', 'secret').start()

        when:
        EnvHelper.withEnv(proxyEnv(proxy, 'foo', 'wrong-pass')) {
            final client = HttpClientHelper
                    .applyProxy(HttpClient.newBuilder().version(HttpClient.Version.HTTP_1_1))
                    .build()
            final request = HttpRequest.newBuilder()
                    .uri(new URI('http://server.proxy-test.internal/hello'))
                    .GET()
                    .build()
            client.send(request, HttpResponse.BodyHandlers.ofString())
        }

        then: 'the client gives up after the proxy rejects the credentials'
        thrown(IOException)
        proxy.unauthorizedCount.get() >= 1
        proxy.proxiedCount.get() == 0

        cleanup:
        proxy?.close()
    }

    @RestoreSystemProperties
    def 'should authenticate a https request tunnelled via proxy connect'() {
        given: 'a self-signed certificate for localhost'
        final keystore = tempDir.resolve('test-keystore.p12')
        def sslContexts = createSelfSignedContexts(keystore)
        and: 'a local https origin server'
        def origin = HttpsServer.create(new InetSocketAddress('127.0.0.1', 0), 0)
        origin.setHttpsConfigurator(new HttpsConfigurator(sslContexts.server))
        origin.createContext('/secure') { exchange ->
            final body = 'secure hello'.getBytes('UTF-8')
            exchange.sendResponseHeaders(200, body.length)
            exchange.getResponseBody().withCloseable { it.write(body) }
        }
        origin.start()
        and: 'an authenticating forward proxy'
        def proxy = new MockAuthProxyServer('foo', 'secret').start()

        when:
        HttpResponse<String> resp = null
        EnvHelper.withEnv(proxyEnv(proxy)) {
            final client = HttpClientHelper
                    .applyProxy(HttpClient.newBuilder()
                        .version(HttpClient.Version.HTTP_1_1)
                        .sslContext(sslContexts.client))
                    .build()
            final request = HttpRequest.newBuilder()
                    .uri(new URI("https://localhost:${origin.address.port}/secure"))
                    .GET()
                    .build()
            resp = client.send(request, HttpResponse.BodyHandlers.ofString())
        }

        then: 'the request is tunnelled through the proxy after the 407 challenge'
        resp.statusCode() == 200
        resp.body() == 'secure hello'
        and:
        proxy.unauthorizedCount.get() >= 1
        proxy.connectCount.get() >= 1

        cleanup:
        proxy?.close()
        origin?.stop(0)
    }

    @RestoreSystemProperties
    def 'should not use the proxy for hosts matching no_proxy'() {
        given: 'a local http origin server'
        def origin = HttpServer.create(new InetSocketAddress('127.0.0.1', 0), 0)
        origin.createContext('/direct') { exchange ->
            final body = 'direct hello'.getBytes('UTF-8')
            exchange.sendResponseHeaders(200, body.length)
            exchange.getResponseBody().withCloseable { it.write(body) }
        }
        origin.start()
        and:
        def proxy = new MockAuthProxyServer('foo', 'secret').start()
        and:
        final env = proxyEnv(proxy)
        env.NO_PROXY = 'example.com,127.0.0.1'

        when:
        HttpResponse<String> resp = null
        EnvHelper.withEnv(env) {
            final client = HttpClientHelper
                    .applyProxy(HttpClient.newBuilder().version(HttpClient.Version.HTTP_1_1))
                    .build()
            final request = HttpRequest.newBuilder()
                    .uri(new URI("http://127.0.0.1:${origin.address.port}/direct"))
                    .GET()
                    .build()
            resp = client.send(request, HttpResponse.BodyHandlers.ofString())
        }

        then: 'the origin is contacted directly, bypassing the proxy'
        resp.statusCode() == 200
        resp.body() == 'direct hello'
        proxy.requestCount.get() == 0

        cleanup:
        proxy?.close()
        origin?.stop(0)
    }

    @RestoreSystemProperties
    def 'should not set any proxy when none is configured'() {
        given:
        System.clearProperty('http.proxyHost')
        System.clearProperty('http.proxyPort')
        System.clearProperty('https.proxyHost')
        System.clearProperty('https.proxyPort')

        when:
        HttpClient client = null
        EnvHelper.withEnv(NO_PROXY_VARS) {
            client = HttpClientHelper.applyProxy(HttpClient.newBuilder()).build()
        }

        then:
        !client.proxy().isPresent()
        !client.authenticator().isPresent()
    }

    private Map<String,SSLContext> createSelfSignedContexts(Path keystore) {
        final password = 'changeit'
        final keytool = Paths.get(System.getProperty('java.home'), 'bin', 'keytool').toString()
        final process = [
                keytool, '-genkeypair',
                '-alias', 'proxy-test',
                '-keyalg', 'RSA', '-keysize', '2048',
                '-validity', '2',
                '-dname', 'CN=localhost',
                '-ext', 'SAN=DNS:localhost,IP:127.0.0.1',
                '-keystore', keystore.toString(),
                '-storetype', 'PKCS12',
                '-storepass', password ].execute()
        assert process.waitFor() == 0, "keytool failed: ${process.err.text}"

        final store = KeyStore.getInstance('PKCS12')
        keystore.withInputStream { store.load(it, password.toCharArray()) }

        final kmf = KeyManagerFactory.getInstance(KeyManagerFactory.getDefaultAlgorithm())
        kmf.init(store, password.toCharArray())
        final serverCtx = SSLContext.getInstance('TLS')
        serverCtx.init(kmf.getKeyManagers(), null, null)

        final tmf = TrustManagerFactory.getInstance(TrustManagerFactory.getDefaultAlgorithm())
        tmf.init(store)
        final clientCtx = SSLContext.getInstance('TLS')
        clientCtx.init(null, tmf.getTrustManagers(), null)

        return [server: serverCtx, client: clientCtx]
    }
}
