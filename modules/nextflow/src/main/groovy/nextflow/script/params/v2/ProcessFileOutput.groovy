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
 * Models a process file output, which defines a file
 * or set of files to be unstaged from a task work directory.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
class ProcessFileOutput {

    /**
     * Lazy expression (e.g. closure) which defines which files
     * to unstage from the task directory.
     * It will be evaluated for each task against the task directory.
     */
    private Object filePattern

    ProcessFileOutput(Object filePattern) {
        this.filePattern = filePattern
    }

    String getFilePattern(Map context) {
        return context.resolveLazy(filePattern)?.toString()
    }

}
