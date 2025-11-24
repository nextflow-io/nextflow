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

package io.seqera.wave.plugin.util

import groovy.transform.Canonical
import groovy.transform.CompileStatic

/**
 * Model a CLI option defined a key-value pairs
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Canonical
@CompileStatic
class BasicCliOpts implements CliOpts {

    Map<String,String> options = new LinkedHashMap<>()
    List<String> args = new ArrayList<>()

    static BasicCliOpts parse(List<String> cli) {
        final result = new BasicCliOpts()
        if( !cli )
            return result
        for( String item : cli ) {
            final p = item.indexOf('=')
            if( p!=-1 ) {
                result.options.put( item.substring(0,p), item.substring(p+1) )
            }
            else {
                result.args.add(item)
            }
        }
        return result
    }
}
