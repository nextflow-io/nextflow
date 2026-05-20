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
package nextflow.processor.hash

import groovy.transform.CompileStatic
import nextflow.SysEnv
import nextflow.processor.TaskRun
/**
 * Factory for creating versioned {@link TaskHasher} instances.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class TaskHasherFactory {

    enum Version {
        STD_V1('std/v1'),
        STD_V2('std/v2')

        final String value

        Version(String value) {
            this.value = value
        }

        String toString() { value }

        static Version of(String val) {
            for( Version v : values() ) {
                if( v.value == val )
                    return v
            }
            throw new IllegalArgumentException("Unknown task hasher version: ${val}")
        }

        static Version DEFAULT() {
            final val = SysEnv.get('NXF_TASK_HASH_VER')
            return val ? of(val) : STD_V2
        }
    }

    static TaskHasher create(TaskRun task) {
        final version = task.processor.session.hashStrategy
        switch( version ) {
            case Version.STD_V1:
                return new TaskHasherV1(task)
            case Version.STD_V2:
                return new TaskHasherV2(task)
            default:
                throw new IllegalArgumentException("Unknown task hasher version: ${version}")
        }
    }
}
