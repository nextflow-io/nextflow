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
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowVariable
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.script.ScriptMeta

/**
 * Exception thrown when trying to invoke a process or function not existing
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class MissingProcessException extends ProcessException {

    MissingProcessException(ScriptMeta meta, MissingMethodException e) {
        super(isChannel(e.type) ? channelMessage(e.type, e.method) : message(meta, e.method, e.arguments), e)
    }

    static boolean isChannel(Class type) {
        return type != null && (DataflowReadChannel.isAssignableFrom(type) || DataflowWriteChannel.isAssignableFrom(type))
    }

    static boolean isValue(Class type) {
        return type != null && DataflowVariable.isAssignableFrom(type)
    }

    static String channelMessage(Class type, String name) {
        return isValue(type)
            ? "Unknown method invocation `$name` on dataflow value -- use an operator such as `map` to transform the value"
            : "Unknown method invocation `$name` on channel type -- use an operator such as `map` to transform channel values"
    }

    @PackageScope
    static String message(ScriptMeta meta, String name, Object[]args) {
        def result = "Missing process or function $name(${args.length ? args : ''})"
        def matches = meta.getAllNames().closest(name)
        if( matches.size()==1 )
            result += " -- Did you mean '${matches[0]}' instead?"
        else if( matches.size()>1 )
            result += "\n\nDid you mean any of these instead?\n" + matches.collect { "  $it"}.join('\n') + '\n'
        return result
    }

}
