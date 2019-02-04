/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

import java.lang.reflect.Field

import groovy.transform.CompileStatic

/**
 * Hack the Java process to be able to access the Unix process id
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class PosixProcess {

    @Delegate
    Process target

    PosixProcess( Process process ) {
        assert process

        target = process
    }

    /**
     * Access the *target* Unix process pid, via a reflection trick
     */
    @Lazy Integer pid = {

        Field field = target.class.getDeclaredField("pid");
        field.setAccessible(true);
        return field.getInt(target);

    } ()


    /**
     * Check if the submitted job has terminated its execution
     */
    boolean hasExited() {

        Field field = target.class.getDeclaredField("hasExited");
        field.setAccessible(true);
        return field.getBoolean(target);

    }

}
