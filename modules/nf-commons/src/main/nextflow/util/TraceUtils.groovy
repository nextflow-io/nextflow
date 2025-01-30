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

import groovy.transform.CompileStatic

/**
 * Implements a bare minimal trace context helper class
 *
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class TraceUtils {

    /**
     * Generates a random {@code traceparent} string to be used for tracing purposes
     *
     * See https://www.w3.org/TR/trace-context/
     *
     * @return A string matching the regexp {@code 00-[a-f0-9]{32}-[a-f0-9]{16}-01}
     */
    static String rndTrace() {
        final traceId = UUID.randomUUID().toString().replace("-", "").substring(0, 32)
        final spanId = UUID.randomUUID().toString().replace("-", "").substring(0, 16)
        return String.format("00-%s-%s-01", traceId, spanId)
    }

}
