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

package nextflow.script.params.v2

import groovy.transform.CompileStatic

/**
 * Models a process file input, which defines a file
 * or set of files to be staged into a task work directory.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
class ProcessFileInput {

    /**
     * File pattern which defines how the input files should be named
     * when they are staged into a task directory.
     */
    private Object filePattern

    /**
     * Lazy expression (e.g. closure) which defines which files
     * to stage in terms of the task inputs.
     * It is evaluated for each task against the task context.
     */
    private Object value

    ProcessFileInput(Object filePattern, Object value) {
        this.filePattern = filePattern != null ? filePattern : '*'
        this.value = value
    }

    String getFilePattern(Map ctx) {
        return ctx.resolveLazy(filePattern)
    }

    Object resolve(Map ctx) {
        return ctx.resolveLazy(value)
    }

}
