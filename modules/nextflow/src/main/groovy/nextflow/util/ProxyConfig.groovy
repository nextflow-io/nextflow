/*
 * Copyright 2013-2024, Seqera Labs
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
 *
 */

package nextflow.util

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode

/**
 * Model a proxy configuration
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * @author Stephen Kazakoff <sh.kazakoff@gmail.com>
 *
 */
@EqualsAndHashCode
@CompileStatic
class ProxyConfig {
    String protocol
    String host
    String port
    String username
    String password

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
                result.username = url.userInfo.substring(0,p)
                result.password = url.userInfo.substring(p+1)
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
}
