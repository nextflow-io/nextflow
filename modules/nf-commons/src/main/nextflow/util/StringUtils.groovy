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

import java.util.regex.Matcher
import java.util.regex.Pattern

import com.google.common.net.InetAddresses
import groovy.transform.CompileStatic

/**
 * String helper routines
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class StringUtils {

    static final public Pattern URL_PROTOCOL = ~/^([a-zA-Z0-9]*):\\/\\/(.+)/
    static final private Pattern URL_PASSWORD = ~/^[a-zA-Z][a-zA-Z0-9]*:\\/\\/(.+)@.+/

    /**
     * Deprecated. Use {@link nextflow.file.FileHelper#getUrlProtocol(java.lang.String)} instead
     */
    @Deprecated
    static String getUrlProtocol(String str) {
        final m = URL_PROTOCOL.matcher(str)
        return m.matches() ? m.group(1) : null
    }

    static final private Pattern BASE_URL = ~/(?i)((?:[a-z][a-zA-Z0-9]*)?:\/\/[^:|\/]+(?::\d*)?)(?:$|\/.*)/

    /**
     * Deprecated. Use {@link nextflow.file.FileHelper#baseUrl(java.lang.String)} instead
     */
    @Deprecated
    static String baseUrl(String url) {
        if( !url )
            return null
        final m = BASE_URL.matcher(url)
        return m.matches() ? m.group(1).toLowerCase() : null
    }

    static private Pattern multilinePattern = ~/["']?(password|token|secret|license)["']?\s?[:=]\s?["']?(\w+)["']?/

    static String stripSecrets(String message) {
        if (message == null) {
            return message
        }
        StringBuilder sb = new StringBuilder(message)
        Matcher matcher = multilinePattern.matcher(sb)
        while (matcher.find()) {
            for(int idx=0; idx<matcher.groupCount(); idx+=2){
                int ini = matcher.start(idx+2)
                int end = matcher.end(idx+2)
                sb.delete(ini, end)
                sb.insert(ini, '********')
            }
            matcher.reset()
        }
        return sb.toString();
    }

    static private boolean isSensitive(Object key) {
        final str = key.toString().toLowerCase()
        return str.contains('password') \
            || str.contains('token') \
            || str.contains('secret') \
            || str.contains('license') \
            || str.contains('auth')
    }


    static Map stripSecrets(Map map) {
        final copy = new HashMap(map.size())
        for( Map.Entry entry : map.entrySet() ) {
            if( entry.value instanceof Map ) {
                copy.put( entry.key, stripSecrets((Map)entry.value))
            }
            else if( isSensitive(entry.key) )  {
                copy.put(entry.key, redact(entry.value))
            }
            else {
                copy.put(entry.key, redactUrlPassword(entry.value))
            }
        }
        return copy
    }

    static String redact(Object value) {
        if( value==null )
            return '(null)'
        if( !value )
            return '(empty)'
        final str = value.toString()
        return str.length()>=10 ? str[0..2] + '****' : '****'
    }

    static String redactUrlPassword(value) {
        final str = value.toString()
        final m = URL_PASSWORD.matcher(str)
        if( m.matches() ) {
            return replaceGroup(m, str, 1, redact(m.group(1)))
        }
        return str
    }

    static String replaceGroup(Matcher matcher, String source, int groupToReplace, String replacement) {
        return new StringBuilder(source)
                .replace(matcher.start(groupToReplace), matcher.end(groupToReplace), replacement)
                .toString()
    }

    static boolean isIpV6String(String address) {
        if( !address || !address.contains(':') )
            return false
        try {
            InetAddresses.forString(address).getAddress().length==16
        }
        catch (IllegalArgumentException e) {
            return false
        }
    }

    static String formatHostName(String host, String port) {
        if( !port || !host )
            return host
        final ipv6 = isIpV6String(host)
        return ipv6 ? "[$host]:$port" : "$host:$port"
    }
}
