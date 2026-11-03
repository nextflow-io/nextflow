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

import io.seqera.http.HxClient
import io.seqera.util.net.ProxyConfig as NetProxyConfig
import spock.lang.Specification
/**
 * Tests for the thin {@link ProxyConfig} adapter that exposes the shared
 * {@code io.seqera:lib-util-net} config to HxClient-based clients. The proxy resolution and
 * matching semantics are covered by {@code io.seqera.util.net.ProxyConfigTest}.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ProxyConfigTest extends Specification {

    def cleanup() {
        ProxyConfig.reset()
    }

    def 'proxyConfig should return null when no proxy is configured'() {
        expect:
        ProxyConfig.proxyConfig() == null
    }

    def 'proxyConfig should return the scheme-aware resolved proxies'() {
        given:
        ProxyConfig.setConfig(NetProxyConfig.fromEnvironment([
                HTTP_PROXY : 'http://http-proxy:3128',
                HTTPS_PROXY: 'http://foo:bar@https-proxy:3129',
                NO_PROXY   : 'internal.example.com' ]))

        when:
        def cfg = ProxyConfig.proxyConfig()

        then:
        cfg != null
        cfg.hasCredentials()
        cfg.toAuthenticator() != null
        and: 'routing is scheme-aware'
        (cfg.toProxySelector().select(URI.create('https://api.example.com'))[0].address() as InetSocketAddress).hostString == 'https-proxy'
        (cfg.toProxySelector().select(URI.create('http://api.example.com'))[0].address() as InetSocketAddress).hostString == 'http-proxy'
        and: 'NO_PROXY entries and loopback bypass the proxy'
        cfg.toProxySelector().select(URI.create('https://internal.example.com')) == [Proxy.NO_PROXY]
        cfg.toProxySelector().select(URI.create('http://localhost:9000')) == [Proxy.NO_PROXY]
    }

    def 'proxyConfig should default the port when missing'() {
        given:
        ProxyConfig.setConfig(NetProxyConfig.fromEnvironment([HTTPS_PROXY: 'https://p']))

        expect:
        (ProxyConfig.proxyConfig().toProxySelector().select(URI.create('https://x/y'))[0].address() as InetSocketAddress).port == 443
    }

    def 'proxyConfig applied via withProxyConfig should reach the client config'() {
        given:
        ProxyConfig.setConfig(NetProxyConfig.fromEnvironment([HTTPS_PROXY: 'http://foo:bar@proxy.example.com:8080']))

        when:
        def client = HxClient.newBuilder().withProxyConfig(ProxyConfig.proxyConfig()).build()

        then:
        client.config.proxySelector != null
        client.config.proxyAuthenticator != null
    }
}
