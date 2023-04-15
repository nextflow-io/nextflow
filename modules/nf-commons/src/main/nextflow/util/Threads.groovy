/*
 * Copyright 2020-2022, Seqera Labs
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
import nextflow.SysEnv
/**
 * Helper class for threads handling
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class Threads {

    static boolean useVirtual() {
        SysEnv.get('NXF_ENABLE_VIRTUAL_THREADS')=='true'
    }

    static Thread start(Closure action) {
        return useVirtual()
                ? Thread.startVirtualThread(action)
                : Thread.startDaemon(action)
    }

    static Thread start(String name, Closure action) {
        if( !useVirtual() )
            return Thread.startDaemon(name, action)
        // create a new virtual thread and change the name
        final result = Thread.startVirtualThread(action)
        result.setName(name)
        return result
    }

}
