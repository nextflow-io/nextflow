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

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString

/**
 * Base CLI options parsed
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@ToString(includeNames = true, includePackage = false)
@EqualsAndHashCode
@CompileStatic
class GnuCliOpts implements CliOpts {
    Map<String,String> options = new LinkedHashMap<>()
    List<String> args = new ArrayList<>()
    String container

    static GnuCliOpts parse(List<String> cli) {
        final result = new GnuCliOpts()

        while( cli.size() ) {
            final name = cli.pop()
            if( name=='--' ) {
                result.args.addAll(cli)
                break
            }
            else if( name.startsWith('-') ) {
                if( cli[0] && !cli[0].startsWith('-') ) {
                    final value = cli.pop()
                    result.options.put(name, value)
                }
                else {
                    result.options.put(name, '')
                }
            }
            else {
                result.container = name
            }
        }
        return result
    }
}
