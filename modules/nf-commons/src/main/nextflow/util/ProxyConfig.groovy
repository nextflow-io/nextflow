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

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import io.seqera.util.net.ProxyConfig as NetProxyConfig

/**
 * Resolve and hold the HTTP/HTTPS/FTP egress proxy for the current run.
 *
 * <p>Resolution and JVM-global installation — the per-protocol {@code <proto>.proxyHost}/{@code .proxyPort}
 * and {@code http.nonProxyHosts} system properties, the default proxy {@link Authenticator} and the
 * Basic-over-{@code CONNECT} tunnelling toggle — is delegated to the shared
 * {@link io.seqera.util.net.ProxyConfig#setupFromEnvironment} in {@code io.seqera:lib-util-net}. The
 * resolved configuration is kept here so HxClient-based clients can be wired via {@link #proxyConfig}.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * @author Stephen Kazakoff <sh.kazakoff@gmail.com>
 */
@Slf4j
@CompileStatic
class ProxyConfig {

    // resolved once by nextflow.cli.Launcher during bootstrap
    private static volatile NetProxyConfig current

    /**
     * Resolve the http/https/ftp proxies from the given environment and install them into the JVM.
     * Called once by {@link nextflow.cli.Launcher} at startup.
     */
    static void setupFromEnvironment(Map<String,String> env) {
        current = NetProxyConfig.setupFromEnvironment(env)
        if( current )
            log.debug "Proxy resolved from environment: $current"
    }

    /**
     * Seed the resolved proxy directly, without the JVM-global installation performed by
     * {@link #setupFromEnvironment} (e.g. to configure it from a source other than the environment).
     */
    static void setConfig(NetProxyConfig cfg) {
        current = cfg
    }

    /**
     * @return The resolved proxy — carrying the proxy selector and (proxy-only) authenticator — or
     *      {@code null} when no proxy is configured. Apply it via
     *      {@link io.seqera.http.HxClient.Builder#withProxyConfig(io.seqera.util.net.ProxyConfig)}.
     */
    static NetProxyConfig proxyConfig() {
        return current
    }

    @PackageScope
    static void reset() {
        current = null
    }
}
