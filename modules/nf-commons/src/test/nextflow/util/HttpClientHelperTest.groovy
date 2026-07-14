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
import java.nio.file.Path
import java.nio.file.Paths
import java.security.KeyStore

import javax.net.ssl.KeyManagerFactory
import javax.net.ssl.SSLContext
import javax.net.ssl.TrustManagerFactory

import com.sun.net.httpserver.HttpsConfigurator
import com.sun.net.httpserver.HttpsServer
import nextflow.SysEnv
import spock.lang.Shared
import spock.lang.Specification
import spock.lang.TempDir
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

    @TempDir
    @Shared
    Path tempDir

    def setup() {
        SysEnv.push(new HashMap<String,String>())
    }

    def cleanup() {
        SysEnv.pop()
    }

    def 'should resolve the proxy configuration from the environment'() {
        when:
        SysEnv.get().putAll(ENV)
        final config = HttpClientHelper.proxyConfigFromEnv()

        then:
        (config != null) == EXPECTED
        config?.hasCredentials() == CREDENTIALS

        where:
        ENV                                                 | EXPECTED | CREDENTIALS
        [:]                                                 | false    | null
        [HTTP_PROXY: 'http://proxy.internal:8080']          | true     | false
        [https_proxy: 'https://proxy.internal:8443']        | true     | false
        [ALL_PROXY: 'http://foo:bar@proxy.internal:8080']   | true     | true
        [HTTPS_PROXY: 'http://foo:bar@proxy.internal:8080'] | true     | true
    }

    def 'should url-decode the proxy credentials'() {
        given:
        SysEnv.get().put('HTTPS_PROXY', 'http://user%40seqera.io:p%40ss%3Aword@proxy.internal:8080')

        when:
        final config = HttpClientHelper.proxyConfigFromEnv()

        then:
        config.hasCredentials()
        and: 'the authenticator releases the decoded credentials for the proxy'
        final pwd = config.toAuthenticator().requestPasswordAuthenticationInstance(
                'proxy.internal', null, 8080, 'http', 'test', 'basic', null,
                Authenticator.RequestorType.PROXY )
        pwd.userName == 'user@seqera.io'
        new String(pwd.password) == 'p@ss:word'
    }

    def 'should authenticate a plain http request against the proxy'() {
        given:
        def proxy = new MockAuthProxyServer('foo', 'secret').start()
        proxy.responseBody = 'hello from proxy'
        SysEnv.get().putAll(proxy.proxyEnv())

        when:
        final resp = sendGet('http://server.proxy-test.internal/hello')

        then: 'the request is answered by the proxy after the 407 challenge'
        resp.statusCode() == 200
        resp.body() == 'hello from proxy'
        and: 'the proxy challenged the client'
        proxy.unauthorizedCount.get() >= 1
        proxy.proxiedCount.get() >= 1

        cleanup:
        proxy?.close()
    }

    def 'should fail when credentials are wrong'() {
        given:
        def proxy = new MockAuthProxyServer('foo', 'secret').start()
        SysEnv.get().putAll(proxy.proxyEnv('foo', 'wrong-pass'))

        when:
        sendGet('http://server.proxy-test.internal/hello')

        then: 'the client gives up after the proxy rejects the credentials'
        thrown(IOException)
        proxy.unauthorizedCount.get() >= 1
        proxy.proxiedCount.get() == 0

        cleanup:
        proxy?.close()
    }

    def 'should authenticate a https request tunnelled via proxy connect'() {
        given: 'a self-signed certificate for the origin host name'
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
        SysEnv.get().putAll(proxy.proxyEnv())

        when: 'the origin is requested via a non-loopback host name, so it is only reachable through the proxy tunnel'
        final resp = sendGet("https://server.proxy-test.internal:${origin.address.port}/secure", sslContexts.client)

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

    def 'should not use the proxy for hosts matching no_proxy'() {
        given:
        SysEnv.get().put('HTTPS_PROXY', 'http://foo:bar@proxy.internal:8080')
        SysEnv.get().put('NO_PROXY', 'example.com')

        when:
        final selector = HttpClientHelper.proxyConfigFromEnv().toProxySelector()

        then: 'hosts matching no_proxy are connected directly'
        selector.select(new URI('https://example.com/foo'))*.type() == [Proxy.Type.DIRECT]
        and: 'loopback addresses are always connected directly'
        selector.select(new URI('http://127.0.0.1/foo'))*.type() == [Proxy.Type.DIRECT]
        selector.select(new URI('https://localhost/foo'))*.type() == [Proxy.Type.DIRECT]
        and: 'any other host goes through the proxy'
        selector.select(new URI('https://seqera.io/foo'))*.type() == [Proxy.Type.HTTP]
    }

    def 'should not set any proxy when none is configured'() {
        when:
        final client = HttpClientHelper.applyProxy(HttpClient.newBuilder()).build()

        then:
        !client.proxy().isPresent()
        !client.authenticator().isPresent()
    }

    /**
     * Send a GET request to the given uri with a client configured via
     * {@link HttpClientHelper#applyProxy}
     */
    private HttpResponse<String> sendGet(String uri, SSLContext sslContext=null) {
        final builder = HttpClient.newBuilder().version(HttpClient.Version.HTTP_1_1)
        if( sslContext )
            builder.sslContext(sslContext)
        final client = HttpClientHelper.applyProxy(builder).build()
        final request = HttpRequest.newBuilder()
                .uri(new URI(uri))
                .GET()
                .build()
        return client.send(request, HttpResponse.BodyHandlers.ofString())
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
                '-ext', 'SAN=DNS:localhost,DNS:server.proxy-test.internal,IP:127.0.0.1',
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
