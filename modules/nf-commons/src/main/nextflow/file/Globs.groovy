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

package nextflow.file

import java.lang.reflect.Method

import groovy.transform.CompileStatic

/**
 * Helper class to access {@link sun.nio.fs.Globs#toRegexPattern(java.lang.String, boolean)} method
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class Globs {

    final private static Method toRegexPattern

    static {
        def clazz = sun.nio.fs.Globs
        toRegexPattern = clazz.getDeclaredMethod('toRegexPattern', String, boolean)
        toRegexPattern.setAccessible(true)
    }

    static String toUnixRegexPattern(String globPattern) {
        toRegexPattern.invoke(null, globPattern, false)
    }

}