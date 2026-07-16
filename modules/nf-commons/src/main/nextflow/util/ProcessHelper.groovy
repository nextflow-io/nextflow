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

import java.lang.management.ManagementFactory

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j

/**
 * Implements helper methods for {@link Process} objects
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class ProcessHelper {

    static long pid(Process process) {
        try {
            // with Java 9 and later pid can be access with `.toHandle().getPid()` method
            final toHandle = Process.class.getMethod("toHandle")
            final handle = toHandle.invoke(process)
            return (Long) Class.forName("java.lang.ProcessHandle").getMethod("pid").invoke(handle)
        }
        catch (Exception e) {
            // ignore it and fallback on next
        }

        // Fallback for Java6+ on unix. This is known not to work on Windows.
        try {
            final pidField = process.class.getDeclaredField("pid")
            pidField.setAccessible(true)
            return pidField.getLong(process)
        }
        catch(Exception e) {
            throw new UnsupportedOperationException("Unable to access process pid for class: ${process.getClass().getName()}",e )
        }
    }

    static long selfPid() {
        try {
            Long.parseLong(ManagementFactory.getRuntimeMXBean().getName().split("@")[0])
        }
        catch(Exception e) {
            throw new UnsupportedOperationException("Unable to find current process pid",e)
        }
    }

}
