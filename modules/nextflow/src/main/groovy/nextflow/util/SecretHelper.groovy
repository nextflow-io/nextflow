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

import java.util.regex.Pattern

import groovy.transform.CompileStatic

/**
 * Helper class to hide sensitive data
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class SecretHelper {

    /**
     * @deprecated Use {@link StringUtils#isSensitiveName(java.lang.String)} instead
     */
    @Deprecated
    static public final Pattern SECRET_KEYS = ~/(?im)^AWS.+|.*TOKEN.*|.*PASSWORD.*|.*SECRET.*|.*accessKey.*|.*api_?key.*/

    /**
     * @deprecated Use {@link #secureEnvString(java.lang.String)} instead
     */
    @Deprecated
    static public final Pattern SECRET_REGEX = ~/(?im)(^AWS[^=]*|.*TOKEN[^=]*|.*SECRET[^=]*|.*API_?KEY[^=]*)=(.*)$/

    // `NAME=` starting only at the beginning of a `[\w.-]` run, possessive, so scanning stays linear
    static private final Pattern ENV_NAME = ~/(?<![\w.-])([\w.-]++)[ \t]*+=/

    /**
     * Mask, up to the end of the line, the value of each `NAME=value` whose name is sensitive
     * by {@link StringUtils#isSensitiveName(java.lang.String)}
     */
    static String secureEnvString( String str ) {
        final result = new StringBuilder()
        final matcher = ENV_NAME.matcher(str)
        int last = 0
        while( matcher.find() ) {
            // a bare `AWS` is not a sensitive config name (it is the `aws` scope) but it was always masked as env var
            final name = matcher.group(1)
            if( !StringUtils.isSensitiveName(name) && !name.equalsIgnoreCase('aws') )
                continue
            int eol = str.indexOf('\n', matcher.end())
            if( eol < 0 )
                eol = str.length()
            result.append(str, last, matcher.end()).append('[secure]')
            last = eol
            matcher.region(eol, str.length())
        }
        return result.append(str, last, str.length()).toString()
    }

    static Object hideSecrets( obj ) {
        if( obj == null )
            return null

        if( obj instanceof Map ) {
            final names = obj.keySet()
            for( String n : names )  {
                // manifest values, e.g. author and license, describe the pipeline and are never credentials
                if( n == 'manifest' )
                    continue
                if( StringUtils.isSensitiveName(n) ) {
                    obj.put(n, '[secret]')
                }
                else {
                    hideSecrets(obj.get(n))
                }
            }
        }
        else if( obj instanceof Collection ) {
            for( Object item : ((Collection)obj) ) {
                hideSecrets(item)
            }
        }
        else if( obj.getClass().isArray() ) {
            for( Object item : ((Object[])obj) ) {
                hideSecrets(item)
            }
        }

        return obj
    }

}
