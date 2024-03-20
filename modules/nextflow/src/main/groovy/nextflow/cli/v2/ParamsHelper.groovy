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
 */

package nextflow.cli.v2

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j

/**
 * Defines the params parsing for commands that accept pipeline params.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class ParamsHelper {

    /**
     * Parse pipeline params from the unmatched args from picocli.
     */
    static Map<String,String> parseParams(List<String> args) {
        Map<String,String> result = [:]
        int i = 0

        while( i < args.size() ) {
            String current = args[i++]
            if( !current.startsWith('--') ) {
                if( current.startsWith('-') )
                    log.warn "Possible legacy command line argument: ${current} -- did you mean -${current} ?"
                throw new IllegalArgumentException("Invalid argument: ${current} -- unable to parse it as a pipeline param or CLI option")
            }

            String key
            String value

            // parse '--param=value'
            if( current.contains('=') ) {
                int split = current.indexOf('=')
                key = current.substring(2, split)
                value = current.substring(split+1)
            }

            // parse '--param value'
            else if( i < args.size() && !args[i].startsWith('--') ) {
                key = current.substring(2)
                value = args[i++]
            }

            // parse '--param1 --param2 ...' as '--param1 true --param2 ...'
            else {
                key = current.substring(2)
                value = 'true'
            }

            result.put(key, value)
        }

        log.trace "Parsing pipeline params from CLI: $result"
        return result
    }

}
