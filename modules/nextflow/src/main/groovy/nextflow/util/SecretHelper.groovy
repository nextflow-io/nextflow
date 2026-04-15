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

    static public final Pattern SECRET_KEYS = ~/(?im)^AWS.+|.*TOKEN.*|.*PASSWORD.*|.*SECRET.*|.*accessKey.*/

    // note: ?i stands for ignore case - ?m stands for multiline
    static public final Pattern SECRET_REGEX = ~/(?im)(^AWS[^=]*|.*TOKEN[^=]*|.*SECRET[^=]*)=(.*)$/

    static String secureEnvString( String str ) {
        str.replaceAll(SECRET_REGEX, '$1=[secure]')
    }

    static Object hideSecrets( obj ) {
        if( obj == null )
            return null

        if( obj instanceof Map ) {
            final names = obj.keySet()
            for( String n : names )  {
                if( SECRET_KEYS.matcher(n).find() ) {
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
