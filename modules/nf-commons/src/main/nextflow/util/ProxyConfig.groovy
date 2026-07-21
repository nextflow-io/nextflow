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

import java.util.concurrent.ConcurrentHashMap

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import io.seqera.http.HxProxyConfig

/**
 * Model a proxy configuration
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * @author Stephen Kazakoff <sh.kazakoff@gmail.com>
 *
 */
@Slf4j
@EqualsAndHashCode
@CompileStatic
class ProxyConfig {
    String protocol
    String host
    String port
    String username
    String password

    // -- proxies resolved from the launch environment, keyed by protocol (http/https/ftp).
    //    Populated once by nextflow.cli.Launcher during bootstrap and exposed via proxyConfig()
    //    so that HxClient-based clients honour the same proxy settings.
    private static final Map<String,ProxyConfig> resolved = new ConcurrentHashMap<>()
    private static volatile List<String> noProxyHosts = List.of()

    @Override
    String toString() {
        def result = protocol ? "protocol=$protocol; host=$host" : "host=$host"
        if( port ) result += "; port=$port"
        if( username ) result += "; username=$username"
        if( password ) result += "; password=${StringUtils.redact(password)}"
        return "ProxyConfig[$result]"
    }

    boolean hasAuthentication() {
        username && password
    }

    Authenticator authenticator() {
        if( !hasAuthentication() )
            return null

        new Authenticator() {
            protected PasswordAuthentication getPasswordAuthentication() {
                if( getRequestorType() != Authenticator.RequestorType.PROXY )
                    return null
                if( !getRequestingHost().equalsIgnoreCase(host) )
                    return null
                if( protocol && !protocol.equalsIgnoreCase(getRequestingProtocol()) )
                    return null
                if( port && getRequestingPort() != Integer.parseInt(port) ) {
                    return null
                }
                return new PasswordAuthentication( username, password.toCharArray() )
            }
        }
    }

    /**
     * Parse a proxy URL string retrieving the host, port, username and password components
     *
     * @param value A proxy string e.g. {@code hostname}, {@code hostname:port}, {@code http://hostname:port},
     *      {@code http://username:password@hostname:port}
     * @return A map object containing at least the host name and, optionally, values for port, username and password.
     *      An empty map if the specified value is empty
     *
     * @throws MalformedURLException when the specified value is not a valid proxy url
     */
    static ProxyConfig parse( final String value ) {
        if( !value )
            return null

        final result = new ProxyConfig()
        int p

        if( value.contains('://') ) {
            def url = new URL(value)
            result.host = url.host
            result.protocol = url.protocol
            if( url.port > 0 )
                result.port = url.port as String
            if( (p=url.userInfo?.indexOf(':') ?: -1) != -1 ) {
                result.username = decodeUserInfo(url.userInfo.substring(0,p))
                result.password = decodeUserInfo(url.userInfo.substring(p+1))
            }
        }
        else if( (p=value.indexOf(':')) != -1 ) {
            result.host = value.substring(0,p)
            result.port = value.substring(p+1)
        }
        else {
            result.host = value
        }

        return result
    }

    /**
     * Register a proxy resolved from the launch environment, so that it can be applied to
     * HxClient-based clients via {@link #proxyConfig}.
     */
    static void register(ProxyConfig proxy) {
        if( proxy?.protocol && proxy.host )
            resolved.put(proxy.protocol.toLowerCase(), proxy)
    }

    @PackageScope
    static void reset() {
        resolved.clear()
        noProxyHosts = List.of()
    }

    /**
     * Record the {@code NO_PROXY} host entries used to bypass the proxy for matching targets.
     */
    static void setNoProxyHosts(List<String> hosts) {
        noProxyHosts = hosts != null
                ? List.copyOf(hosts.collect { it.trim().toLowerCase() }.findAll { it } as List<String>)
                : List.of()
    }

    /**
     * @return An {@link HxProxyConfig} view of the resolved proxy environment — the HxClient
     *      library type carrying the proxy selector and (proxy-only) authenticator — or
     *      {@code null} when no proxy is configured. Apply it via
     *      {@link HxClient.Builder#withProxyConfig(HxProxyConfig)}.
     */
    static HxProxyConfig proxyConfig() {
        final http = resolved.get('http')
        final https = resolved.get('https')
        if( !http && !https )
            return null
        final builder = HxProxyConfig.newBuilder()
        if( http )
            builder.httpProxy(http.host, portAsInt(http.port, 80), http.username, http.password)
        if( https )
            builder.httpsProxy(https.host, portAsInt(https.port, 443), https.username, https.password)
        builder.noProxy(noProxyHosts)
        return builder.build()
    }

    /**
     * Percent-decode a userinfo component (username or password) per RFC 3986, so proxy
     * credentials carrying special characters (e.g. {@code @}, {@code :}) work. A literal
     * {@code +} is preserved — userinfo is not form-encoded — so it is shielded from the
     * {@code +}→space rule of {@link URLDecoder}.
     */
    private static String decodeUserInfo(String s) {
        return s != null ? URLDecoder.decode(s.replace('+', '%2B'), 'UTF-8') : null
    }

    private static int portAsInt(String port, int defaultPort) {
        if( !port )
            return defaultPort
        try {
            return Integer.parseInt(port.trim())
        }
        catch( NumberFormatException e ) {
            log.warn("Ignoring invalid proxy port '$port' - using default $defaultPort")
            return defaultPort
        }
    }
}
