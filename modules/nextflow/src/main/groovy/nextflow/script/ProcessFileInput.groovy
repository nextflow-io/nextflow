/*
 * Copyright 2013-2023, Seqera Labs
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

package nextflow.script

import groovy.transform.CompileStatic
import nextflow.util.LazyHelper

/**
 * Models a process file input, which defines a file
 * or set of files to be staged into a task work directory.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
class ProcessFileInput implements PathArityAware {

    private Object value

    private String name

    /**
     * Flag to support legacy `file` input.
     */
    private boolean pathQualifier

    private Object filePattern

    ProcessFileInput(Object value, String name, boolean pathQualifier, Map<String,?> opts) {
        this.value = value
        this.name = name
        this.pathQualifier = pathQualifier
        this.filePattern = opts.stageAs ?: opts.name

        for( Map.Entry<String,?> entry : opts )
            setProperty(entry.key, entry.value)
    }

    void setStageAs(String value) {
        this.filePattern = value
    }

    Object resolve(Map ctx) {
        return LazyHelper.resolve(ctx, value)
    }

    String getName() {
        return name
    }

    boolean isPathQualifier() {
        return pathQualifier
    }

    String getFilePattern(Map ctx) {
        if( filePattern != null )
            return LazyHelper.resolve(ctx, filePattern)
        else
            return filePattern = '*'
    }

}
