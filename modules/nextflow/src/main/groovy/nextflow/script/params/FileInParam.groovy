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

package nextflow.script.params

import groovy.transform.InheritConstructors
import groovy.util.logging.Slf4j

/**
 * Represents a process *file* input parameter
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@InheritConstructors
class FileInParam extends BaseInParam implements PathQualifier {

    private boolean pathQualifier

    private Map<String,?> options

    @Override String getTypeName() { pathQualifier ? 'path' : 'file' }

    @Override String getTypeSimpleName() { getTypeName() + "inparam" }

    String getName() {
        if( bindObject instanceof Map ) {
            assert !pathQualifier
            def entry = bindObject.entrySet().first()
            return entry?.key
        }

        if( bindObject instanceof GString ) {
            return '__$' + this.toString()
        }

        return super.getName()
    }

    @Override
    BaseInParam bind( obj ) {
        if( pathQualifier && obj instanceof Map )
            throw new IllegalArgumentException("Input `path` does not allow such arguments: ${obj.entrySet().collect{"${it.key}:${it.value}"}.join(',')}")
        super.bind(obj)
        return this
    }

    @Override
    FileInParam setPathQualifier(boolean flag) {
        pathQualifier = flag
        return this
    }

    @Override
    boolean isPathQualifier() { pathQualifier }

    @Override
    FileInParam setOptions(Map<String,?> opts) {
        this.options = opts
        return this
    }

    Map<String,?> getOptions() { options }

}
