/*
 * Copyright 2020-2021, Seqera Labs
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

import java.util.regex.Pattern

import groovy.transform.CompileStatic

/**
 * String helper routines
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class StringUtils {

    static final public Pattern URL_PROTOCOL = ~/^([a-zA-Z0-9]*):\\/\\/(.+)/

    static String getUrlProtocol(String str) {
        final m = URL_PROTOCOL.matcher(str)
        return m.matches() ? m.group(1) : null
    }

    static final private Pattern BASE_URL = ~/(?i)((?:[a-z][a-zA-Z0-9]*)?:\/\/[^:|\/]+(?::\d*)?)(?:$|\/.*)/

    static String baseUrl(String url) {
        if( !url )
            return null
        final m = BASE_URL.matcher(url)
        return m.matches() ? m.group(1).toLowerCase() : null
    }
}
