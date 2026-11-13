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

package nextflow.exception

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import nextflow.script.ScriptMeta

/**
 * Exception thrown when a module component cannot be found
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class MissingModuleComponentException extends ProcessException {

    MissingModuleComponentException(ScriptMeta meta, String name) {
        super(message(meta, name))
    }

    @PackageScope
    static String message(ScriptMeta meta, String name) {
        def result = "Cannot find a component with name '$name' in module: $meta.scriptPath"
        def names = meta.getDefinitions().findAll { it.name }.collect { it.name }
        def matches = names.closest(name)
        if( matches )
            result += "\n\nDid you mean any of these?\n" + matches.collect { "  $it"}.join('\n') + '\n'
        return result
    }
}
