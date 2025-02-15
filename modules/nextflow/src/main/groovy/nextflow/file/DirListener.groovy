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

package nextflow.file

import static java.nio.file.StandardWatchEventKinds.*

import java.nio.file.Path
import java.nio.file.WatchEvent

import groovy.transform.stc.ClosureParams
import groovy.transform.stc.SimpleType
/**
 * Declare interface methods for directory watching events
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
interface DirListener {

    static public Map<String, WatchEvent.Kind<Path>> EVENT_MAP = [
            'create':ENTRY_CREATE,
            'delete':ENTRY_DELETE,
            'modify':ENTRY_MODIFY
    ]

    void apply(@ClosureParams(value = SimpleType, options = ['java.nio.file.Path']) Closure onNext )

    void onComplete(Closure action)

    void terminate()

    /**
     * Converts a comma separated events string to the corresponding {@code WatchEvent.Kind} instances
     *
     * @param events the list of events to watch
     * @return
     */
    default WatchEvent.Kind<Path>[] stringToWatchEvents(String events = null){
        def result = []
        if( !events )
            result << ENTRY_CREATE

        else {
            events.split(',').each { String it ->
                def ev = it.trim().toLowerCase()
                def val = EVENT_MAP[ev]
                if( !val )
                    throw new IllegalArgumentException("Invalid watch event: $it -- Valid values are: ${EVENT_MAP.keySet().join(', ')}")
                result << val
            }
        }

        result as WatchEvent.Kind<Path>[]
    }
}
