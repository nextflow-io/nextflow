/*
 * Copyright 2013-2023, Seqera Labs
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
import groovy.transform.PackageScope
import org.apache.commons.lang.StringUtils

/**
 * A {@link Map} that automatically converts kebab-case
 * keys to camelCase.
 * 
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
class AliasMap implements Map<String,Object> {

    @Delegate
    private Map<String,Object> target

    AliasMap() {
        target = new LinkedHashMap<>()
    }

    AliasMap(Map<String,Object> values) {
        this()
        putAll(values)
    }

    @Override
    Object put(String key, Object value) {
        if( key.contains('-') )
            key = kebabToCamelCase(key)

        return target.put(key, value)
    }

    @PackageScope
    static String kebabToCamelCase( String str ) {
        if( !str )
            return str

        final result = new StringBuilder()
        str.split('-').eachWithIndex { String entry, int i ->
            result << (i>0 ? StringUtils.capitalize(entry) : entry )
        }

        return result.toString()
    }

    @PackageScope
    static String camelToKebabCase( String str ) {
        final lower = 'a'..'z'
        final upper = 'A'..'Z'

        final result = new StringBuilder()
        if( !str )
            return str

        result << str[0]
        for( int i=1; i<str.size(); i++ ) {
            if( str[i] in upper && str[i-1] in lower  ) {
                result << '-'
                if( i+1<str.size() && str[i+1] in lower ) {
                    result << str[i].toLowerCase()
                }
                else {
                    result << str[i]
                }
            }
            else {
                result << str[i]
            }
        }

        return result.toString()
    }

}
