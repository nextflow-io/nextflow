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

package nextflow.util

/**
 * Some utils for debugging
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
final class DebugUtils {

    private static final String SEPARATOR = "\n";

    private DebugUtils() {
    }

    static String formatStackTrace(StackTraceElement[] stackTrace) {
        StringBuilder buffer = new StringBuilder();
        for (StackTraceElement element : stackTrace) {
            buffer.append(element).append(SEPARATOR);
        }
        return buffer.toString();
    }

    static String formatCurrentStacktrace() {
        StackTraceElement[] stackTrace = Thread.currentThread().getStackTrace();
        return formatStackTrace(stackTrace);
    }
}
