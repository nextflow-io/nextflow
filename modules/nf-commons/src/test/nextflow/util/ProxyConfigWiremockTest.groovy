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

import java.net.http.HttpRequest
import java.net.http.HttpResponse

import com.github.tomakehurst.wiremock.WireMockServer
import com.github.tomakehurst.wiremock.client.WireMock
import io.seqera.http.HxClient
import spock.lang.Specification

import static com.github.tomakehurst.wiremock.client.WireMock.*
import static com.github.tomakehurst.wiremock.core.WireMockConfiguration.wireMockConfig

/**
 * Integration test proving that an {@link HxClient} configured from {@link ProxyConfig} actually
 * routes a request through the forward proxy and carries the proxy credentials.
 *
 * <p>WireMock stands in for the proxy: for a plain-HTTP target the JDK client issues the request in
 * absolute-URI form directly to the proxy, so WireMock receives it and can enforce the
 * {@code Proxy-Authorization} challenge/response handshake (no HTTPS CONNECT tunnelling involved,
 * hence no {@code jdk.http.auth.tunneling.disabledSchemes} caveat).
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ProxyConfigWiremockTest extends Specification {

    WireMockServer wireMock

    def setup() {
        wireMock = new WireMockServer(wireMockConfig().dynamicPort())
        wireMock.start()
        WireMock.configureFor('localhost', wireMock.port())
    }

    def cleanup() {
        wireMock?.stop()
        ProxyConfig.reset()
    }

    private void stubAuthenticatingProxy() {
        // authorized request -> 200
        stubFor(get(urlPathEqualTo('/hello')).atPriority(1)
                .withHeader('Proxy-Authorization', matching('Basic .*'))
                .willReturn(aResponse().withStatus(200).withBody('OK via proxy')))
        // unauthorized request -> 407 challenge
        stubFor(get(urlPathEqualTo('/hello')).atPriority(2)
                .willReturn(aResponse().withStatus(407).withHeader('Proxy-Authenticate', 'Basic realm="test"')))
    }

    private static HttpRequest request() {
        // the target is a non-loopback host so it is not bypassed; with a proxy configured the JDK
        // never resolves it — the request is sent to the proxy in absolute-URI form
        return HttpRequest.newBuilder().uri(URI.create('http://proxy-target.example/hello')).GET().build()
    }

    def 'should route an HxClient request through the authenticating proxy'() {
        given: 'a resolved http proxy with credentials pointing at WireMock'
        stubAuthenticatingProxy()
        ProxyConfig.register(new ProxyConfig(protocol: 'http', host: '127.0.0.1', port: wireMock.port() as String, username: 'alice', password: 's3cret'))

        and: 'an HxClient wired the same way the production sites are'
        def client = HxClient.newBuilder()
                .withProxyConfig(ProxyConfig.proxyConfig())
                .build()

        when:
        def response = client.send(request(), HttpResponse.BodyHandlers.ofString())

        then: 'the request succeeded through the proxy after the credential handshake'
        response.statusCode() == 200
        response.body() == 'OK via proxy'

        and: 'the proxy received a request carrying the Basic credentials'
        wireMock.verify(getRequestedFor(urlPathEqualTo('/hello')).withHeader('Proxy-Authorization', matching('Basic .*')))
    }

    def 'should get 407 when the proxy requires credentials that are not configured'() {
        given: 'a resolved http proxy WITHOUT credentials'
        stubAuthenticatingProxy()
        ProxyConfig.register(new ProxyConfig(protocol: 'http', host: '127.0.0.1', port: wireMock.port() as String))

        and:
        def client = HxClient.newBuilder()
                .withProxyConfig(ProxyConfig.proxyConfig())
                .maxAttempts(1)
                .build()

        when:
        def response = client.send(request(), HttpResponse.BodyHandlers.ofString())

        then: 'routing worked (the request reached the proxy) but the challenge was unanswered'
        response.statusCode() == 407
        wireMock.verify(getRequestedFor(urlPathEqualTo('/hello')))
    }
}
