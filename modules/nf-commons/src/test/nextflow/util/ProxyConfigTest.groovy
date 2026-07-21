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
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ProxyConfigTest extends Specification {

    def cleanup() {
        ProxyConfig.reset()
    }

    private static PasswordAuthentication invokeAuth(Authenticator authenticator, Authenticator.RequestorType type, String host, int port, String protocol) {
        // drive the protected Authenticator API the same way the JDK HTTP client does
        return authenticator.requestPasswordAuthenticationInstance(
                host, null, port, protocol, 'proxy auth', 'basic', new URL('http://example.com'), type )
    }

    def 'should parse proxy env variables'( ) {

        expect:
        ProxyConfig.parse(null) == null

        ProxyConfig.parse('http://domain') == new ProxyConfig(protocol: 'http', host: 'domain')
        ProxyConfig.parse('http://domain:333') == new ProxyConfig(protocol: 'http',host: 'domain', port: '333')
        ProxyConfig.parse('http://10.20.30.40') == new ProxyConfig(protocol: 'http',host: '10.20.30.40')
        ProxyConfig.parse('http://10.20.30.40:333') == new ProxyConfig(protocol: 'http',host: '10.20.30.40', port: '333')
        ProxyConfig.parse('http://10.20.30.40:333/some/path') == new ProxyConfig(protocol: 'http',host: '10.20.30.40', port: '333')

        ProxyConfig.parse('http://user:pass@domain') == new ProxyConfig(protocol: 'http',host: 'domain', username: 'user', password: 'pass')
        ProxyConfig.parse('http://user:pass@domain:333') == new ProxyConfig(protocol: 'http',host: 'domain', port: '333', username: 'user', password: 'pass')
        ProxyConfig.parse('http://user:pass@10.20.30.40') == new ProxyConfig(protocol: 'http',host: '10.20.30.40', username: 'user', password: 'pass')
        ProxyConfig.parse('http://user:pass@10.20.30.40:333') == new ProxyConfig(protocol: 'http',host: '10.20.30.40', port: '333', username: 'user', password: 'pass')
        ProxyConfig.parse('http://user:pass@10.20.30.40:333/some/path') == new ProxyConfig(protocol: 'http',host: '10.20.30.40', port: '333', username: 'user', password: 'pass')

        ProxyConfig.parse('foo') == new ProxyConfig(host: 'foo')
        ProxyConfig.parse('foo:123') == new ProxyConfig(host: 'foo', port: '123')

    }

    def 'should url-decode proxy credentials'() {
        expect: 'percent-encoded username/password are decoded'
        ProxyConfig.parse('http://user%40corp:p%40ss%3Aword@domain:333') == new ProxyConfig(protocol: 'http', host: 'domain', port: '333', username: 'user@corp', password: 'p@ss:word')
        and: 'a plus sign stays literal (RFC 3986 userinfo, not form encoding)'
        ProxyConfig.parse('http://a+b:c+d@domain') == new ProxyConfig(protocol: 'http', host: 'domain', username: 'a+b', password: 'c+d')
    }

    def 'authenticator should release credentials only for matching proxy challenges'() {
        given:
        def proxy = new ProxyConfig(protocol: 'https', host: 'proxy.example.com', port: '8080', username: 'foo', password: 'bar')
        def authenticator = proxy.authenticator()

        expect: 'proxy challenge from the matching host/port/protocol'
        with(invokeAuth(authenticator, Authenticator.RequestorType.PROXY, 'proxy.example.com', 8080, 'https')) {
            userName == 'foo'
            new String(password) == 'bar'
        }
        and: 'origin-server challenge yields nothing'
        invokeAuth(authenticator, Authenticator.RequestorType.SERVER, 'proxy.example.com', 8080, 'https') == null
        and: 'different host yields nothing'
        invokeAuth(authenticator, Authenticator.RequestorType.PROXY, 'other.example.com', 8080, 'https') == null
    }

    def 'authenticator should release credentials for an https proxy over an http CONNECT tunnel'() {
        given: 'an https_proxy-derived config (Launcher sets protocol=https from HTTPS_PROXY)'
        def proxy = new ProxyConfig(protocol: 'https', host: 'proxy.example.com', port: '8080', username: 'foo', password: 'bar')
        def authenticator = proxy.authenticator()

        when: 'the JDK authenticates an HTTPS CONNECT tunnel — it reports requestingProtocol="http", not "https" (see #5634)'
        def auth = invokeAuth(authenticator, Authenticator.RequestorType.PROXY, 'proxy.example.com', 8080, 'http')

        then: 'credentials must still be released, otherwise HTTPS targets 407 behind an authenticating proxy'
        auth != null
        auth.userName == 'foo'
        new String(auth.password) == 'bar'
    }

    def 'authenticator relaxation is one-directional: http proxy config rejects an https request'() {
        given: 'an http_proxy-derived config'
        def proxy = new ProxyConfig(protocol: 'http', host: 'proxy.example.com', port: '8080', username: 'foo', password: 'bar')
        def authenticator = proxy.authenticator()

        expect: 'an https requesting-protocol against an http config is still rejected'
        invokeAuth(authenticator, Authenticator.RequestorType.PROXY, 'proxy.example.com', 8080, 'https') == null
    }

    def 'authenticator is null-safe when the requesting protocol is missing'() {
        given: 'an https_proxy-derived config'
        def proxy = new ProxyConfig(protocol: 'https', host: 'proxy.example.com', port: '8080', username: 'foo', password: 'bar')
        def authenticator = proxy.authenticator()

        expect: 'a null requesting protocol does not throw and yields no credentials'
        invokeAuth(authenticator, Authenticator.RequestorType.PROXY, 'proxy.example.com', 8080, null) == null
    }

    def 'proxyConfig should return null when no proxy is registered'() {
        expect:
        ProxyConfig.proxyConfig() == null
    }

    def 'proxyConfig should build an HxProxyConfig from the resolved proxies'() {
        given:
        ProxyConfig.register(new ProxyConfig(protocol: 'http', host: 'http-proxy', port: '3128'))
        ProxyConfig.register(new ProxyConfig(protocol: 'https', host: 'https-proxy', port: '3129', username: 'foo', password: 'bar'))
        ProxyConfig.setNoProxyHosts(['internal.example.com'])

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
        ProxyConfig.register(new ProxyConfig(protocol: 'https', host: 'p'))

        expect:
        (ProxyConfig.proxyConfig().toProxySelector().select(URI.create('https://x/y'))[0].address() as InetSocketAddress).port == 443
    }

    def 'proxyConfig applied via withProxyConfig should reach the client config'() {
        given:
        ProxyConfig.register(new ProxyConfig(protocol: 'https', host: 'proxy.example.com', port: '8080', username: 'foo', password: 'bar'))

        when:
        def client = HxClient.newBuilder().withProxyConfig(ProxyConfig.proxyConfig()).build()

        then:
        client.config.proxySelector != null
        client.config.proxyAuthenticator != null
    }
}
