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

package io.seqera.wave.plugin

import java.net.http.HttpRequest
import java.net.http.HttpResponse

import nextflow.Session
import nextflow.SysEnv
import spock.lang.Specification
import test.MockAuthProxyServer

/**
 * Verify the Wave client HTTP transport authenticates against a forward proxy
 * requiring Basic authentication resolved from the environment.
 *
 * @author Phil Ewels <phil.ewels@seqera.io>
 */
class WaveClientProxyTest extends Specification {

    def 'should send wave requests via an authenticated proxy'() {
        given: 'an authenticating forward proxy'
        def proxy = new MockAuthProxyServer('foo', 'secret').start()
        proxy.responseContentType = 'application/json'
        proxy.responseBody = '{"status":"OK"}'
        and: 'the proxy settings with credentials defined in the environment'
        SysEnv.push(proxy.proxyEnv())

        when:
        final session = Mock(Session) { getConfig() >> [:] }
        final wave = new WaveClient(session)
        final client = wave.newHttpClient0()
        final request = HttpRequest.newBuilder()
                .uri(new URI('http://wave.proxy-test.internal/service-info'))
                .GET()
                .build()
        final resp = client.send(request, HttpResponse.BodyHandlers.ofString())

        then: 'the client is configured with the proxy and its credentials'
        client.proxy().isPresent()
        client.authenticator().isPresent()
        and: 'the request is answered by the proxy after the 407 challenge'
        resp.statusCode() == 200
        resp.body() == '{"status":"OK"}'
        and:
        proxy.unauthorizedCount.get() >= 1
        proxy.proxiedCount.get() >= 1

        cleanup:
        proxy?.close()
        SysEnv.pop()
    }
}
