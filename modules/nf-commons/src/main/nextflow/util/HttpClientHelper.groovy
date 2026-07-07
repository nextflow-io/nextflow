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
import io.seqera.http.HxProxyConfig

/**
 * Helper methods for {@link HttpClient} instances
 *
 * @author Phil Ewels <phil.ewels@seqera.io>
 */
@Slf4j
@CompileStatic
class HttpClientHelper {

    /**
     * Apply the forward proxy settings resolved from the environment to the
     * given {@link HttpClient.Builder}.
     *
     * The proxy is resolved from the {@code HTTPS_PROXY}, {@code HTTP_PROXY}
     * and {@code ALL_PROXY} environment variables (both cases) with Java system
     * properties as fallback, supporting URL-encoded {@code user:pass@host:port}
     * credentials. Hosts matching {@code NO_PROXY} are connected directly.
     *
     * Note: this is required because {@link HttpClient} ignores the JVM default
     * authenticator ({@link java.net.Authenticator#setDefault}) — proxy credentials
     * take effect only when the authenticator is set on the client builder itself.
     * Clients built via a plain {@code HxClient.newBuilder()} resolve the proxy
     * automatically; this helper is only needed for pre-built {@link HttpClient}
     * instances passed via {@code HxClient.Builder.httpClient(...)}.
     *
     * @param builder The {@link HttpClient.Builder} to configure
     * @return The same builder instance, with the proxy selector and proxy
     *      authenticator applied when a proxy is configured in the environment
     */
    static HttpClient.Builder applyProxy(HttpClient.Builder builder) {
        final proxyConfig = HxProxyConfig.fromEnvironment()
        if( proxyConfig ) {
            log.trace "Applying proxy settings to HTTP client: host=${proxyConfig.host()}; port=${proxyConfig.port()}; credentials=${proxyConfig.hasCredentials()}"
            builder.proxy(proxyConfig.toProxySelector())
            final authenticator = proxyConfig.toAuthenticator()
            if( authenticator )
                builder.authenticator(authenticator)
        }
        return builder
    }
}
