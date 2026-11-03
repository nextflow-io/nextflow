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

    // the `name:` or `name=` of a pair; the name starts only at the beginning of a `[\w-]` run and
    // possessive quantifiers never backtrack, so scanning stays linear on long inputs
    static private final Pattern KEY_SEP = ~/(?<![\w-])([\w-]++)["']?\s*+[:=]\s*+/

    // the value of a sensitive pair; an unquoted value ends at a quote, a space, `]` or `, `
    // (the entry separator of a map rendered as text)
    static private final Pattern VALUE = ~/(?i:bearer|basic)[ \t]++(?:[^\s"',\]]|,(?!\s))++|"[^"]*+"|'[^']*+'|["']?(?:[^\s"',\]]|,(?!\s))++/

    // stricter than isSensitiveName: in free text any word followed by `:` or `=` is a candidate
    static private final Pattern SENSITIVE_WORD = ~/(?i)auth\d*$|authorization|token|password|passwd|secret|license|credential|(?:_key|accountkey|accesskey|apikey|privatekey)$/

    static private boolean isSensitiveWord(String word) {
        return !word.equalsIgnoreCase('tokens') && SENSITIVE_WORD.matcher(word).find()
    }

    // lower case parts of names that may hold a credential; `apikey` covers LLM provider keys, e.g. OPENAI_API_KEY
    static public final List<String> SENSITIVE_NAME_PARTS = ['password','passwd','token','secret','license','auth','credential',
        'apikey','api_key','accesskey','access_key','privatekey','private_key','accountkey','account_key'].asImmutable()

    /**
     * @return {@code true} when a config option or environment variable name may hold a credential
     */
    static boolean isSensitiveName(String name) {
        if( !name )
            return false
        final str = name.toLowerCase()
        return (str.length()>3 && str.startsWith('aws')) || str.endsWith('_key') || SENSITIVE_NAME_PARTS.any { str.contains(it) }
    }

    static String stripSecrets(String message) {
        if( message == null )
            return message
        // the value is read only for a sensitive name, so a non-sensitive value never hides a pair that follows it
        final result = new StringBuilder()
        final key = KEY_SEP.matcher(message).useTransparentBounds(true)
        final value = VALUE.matcher(message)
        int last = 0
        while( key.find() ) {
            if( !isSensitiveWord(key.group(1).replace('-','_')) )
                continue
            value.region(key.end(), message.length())
            if( !value.lookingAt() )
                continue
            final str = value.group()
            final open = str.startsWith('"') || str.startsWith("'") ? str.substring(0,1) : ''
            final close = open && str.length()>1 && str.endsWith(open) ? open : ''
            result.append(message, last, key.end()).append(open).append('********').append(close)
            last = value.end()
            key.region(last, message.length())
        }
        return result.append(message, last, message.length()).toString()
    }


    static Map stripSecrets(Map map) {
        if( map==null )
            return null
        final copy = new HashMap(map.size())
        for( Map.Entry entry : map.entrySet() ) {
            if( entry.value instanceof Map ) {
                copy.put( entry.key, stripSecrets((Map)entry.value))
            }
            else if( isSensitiveName(entry.key?.toString()) )  {
                copy.put(entry.key, entry.value ? '****' : redact(entry.value))
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
