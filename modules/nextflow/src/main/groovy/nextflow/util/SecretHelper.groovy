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

    // NOTE: `api_?key` covers `agent.apiKey` and the `*_API_KEY` environment spellings; without it
    // a resolved provider credential is persisted verbatim by the lineage observer and shipped as
    // `workflow.configText` (`accessKey` does NOT match `apiKey`)
    static public final Pattern SECRET_KEYS = ~/(?im)^AWS.+|.*TOKEN.*|.*PASSWORD.*|.*SECRET.*|.*accessKey.*|.*api_?key.*/

    // note: ?i stands for ignore case - ?m stands for multiline
    // NOTE: `api_?key` mirrors SECRET_KEYS above. The two must agree: SECRET_KEYS masks a config
    // MAP entry while this masks a `NAME=value` environment line, and the out-of-band credential
    // channel the agent docs still recommend (`env { OPENAI_API_KEY = ... }`,
    // `agent.containerOptions = '-e OPENAI_API_KEY'`) is delivered as exactly such a line -- so a
    // pattern that covers only one of the two masks the credential in only half the places it
    // appears. Kept in sync with the identical twin in {@code nextflow.trace.TraceRecord}.
    static public final Pattern SECRET_REGEX = ~/(?im)(^AWS[^=]*|.*TOKEN[^=]*|.*SECRET[^=]*|.*API_?KEY[^=]*)=(.*)$/

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
