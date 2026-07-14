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

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import io.seqera.http.HxClient
import io.seqera.http.HxProxyConfig
import nextflow.SysEnv

/**
 * Helper methods to configure HTTP clients with the forward proxy settings
 * resolved from the environment.
 *
 * The proxy is resolved from the {@code HTTPS_PROXY}, {@code HTTP_PROXY} and
 * {@code ALL_PROXY} environment variables (upper and lower case), supporting
 * URL-encoded {@code user:pass@host:port} credentials. Hosts matching
 * {@code NO_PROXY} are connected directly. Resolving the environment here is
 * by design: lib-httpx expects the proxy configuration to be provided
 * explicitly via its API (see seqeralabs/libseqera#82).
 *
 * Note: this is required because {@link HttpClient} ignores the JVM default
 * authenticator ({@link java.net.Authenticator#setDefault}) — proxy credentials
 * take effect only when the authenticator is set on the client builder itself.
 *
 * @author Phil Ewels <phil.ewels@seqera.io>
 */
@Slf4j
@CompileStatic
class HttpClientHelper {

    /**
     * Resolve the forward proxy configuration from the environment.
     *
     * @return A {@link HxProxyConfig} holding the proxy settings, or {@code null}
     *      when no proxy is configured in the environment
     */
    static HxProxyConfig proxyConfigFromEnv() {
        final env = SysEnv.get()
        final allProxy = env.get('ALL_PROXY') ?: env.get('all_proxy')
        final httpProxy = parse(env.get('HTTP_PROXY') ?: env.get('http_proxy') ?: allProxy)
        final httpsProxy = parse(env.get('HTTPS_PROXY') ?: env.get('https_proxy') ?: allProxy)
        if( !httpProxy && !httpsProxy )
            return null
        final builder = HxProxyConfig.newBuilder()
        if( httpProxy )
            builder.httpProxy(httpProxy.host, portOf(httpProxy), httpProxy.username, httpProxy.password)
        if( httpsProxy )
            builder.httpsProxy(httpsProxy.host, portOf(httpsProxy), httpsProxy.username, httpsProxy.password)
        final noProxy = env.get('NO_PROXY') ?: env.get('no_proxy')
        if( noProxy )
            builder.noProxy(noProxy.tokenize(',')*.trim())
        return builder.build()
    }

    private static ProxyConfig parse(String value) {
        try {
            return ProxyConfig.parse(value)
        }
        catch( MalformedURLException e ) {
            log.warn "Not a valid proxy URL: '$value' -- Check the value of the proxy variables in your environment"
            return null
        }
    }

    private static int portOf(ProxyConfig proxy) {
        if( proxy.port )
            return proxy.port as int
        return proxy.protocol == 'https' ? 443 : 80
    }

    /**
     * Apply the forward proxy settings resolved from the environment to the
     * given {@link HttpClient.Builder}. When no proxy is configured the builder
     * is left untouched, preserving the JDK default proxy selector behaviour.
     *
     * @param builder The {@link HttpClient.Builder} to configure
     * @return The same builder instance, with the proxy selector and proxy
     *      authenticator applied when a proxy is configured in the environment
     */
    static HttpClient.Builder applyProxy(HttpClient.Builder builder) {
        final config = proxyConfigFromEnv()
        if( config ) {
            log.trace "Applying proxy settings to HTTP client: credentials=${config.hasCredentials()}"
            builder.proxy(config.toProxySelector())
            final authenticator = config.toAuthenticator()
            if( authenticator )
                builder.authenticator(authenticator)
        }
        return builder
    }

    /**
     * Apply the forward proxy settings resolved from the environment to the
     * given {@link HxClient.Builder}. When no proxy is configured the builder
     * is left untouched, preserving the JDK default proxy selector behaviour.
     *
     * @param builder The {@link HxClient.Builder} to configure
     * @return The same builder instance, with the proxy configuration applied
     *      when a proxy is configured in the environment
     */
    static HxClient.Builder applyProxy(HxClient.Builder builder) {
        final config = proxyConfigFromEnv()
        if( config ) {
            log.trace "Applying proxy settings to HTTP client: credentials=${config.hasCredentials()}"
            builder.withProxyConfig(config)
        }
        return builder
    }
}
