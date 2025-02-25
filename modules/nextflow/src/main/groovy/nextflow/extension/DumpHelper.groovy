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

package nextflow.extension

import java.nio.file.Path

import groovy.json.JsonOutput
import groovy.transform.CompileStatic
import nextflow.extension.FilesEx
import org.codehaus.groovy.runtime.InvokerHelper
/**
 * Helper converters for {@link DumpOp} operator
 *
 * @author : jorge <jorge.aguilera@seqera.io>
 */
@CompileStatic
class DumpHelper {

    static def deepConvertToString(value, nullValue = '') {
        if( value instanceof List )
            value.collect { it -> deepConvertToString(it) }

        else if( value instanceof Map )
            value.inject([:]) { accum, it ->
                accum[it.key] = deepConvertToString(it.value)
                accum
            }

        else if( value instanceof Path )
            FilesEx.toUriString((Path)value)

        else if( value == null )
            nullValue

        else
            value.toString()
    }

    static String prettyPrint(value) {
        if( value instanceof List || value instanceof Map || value instanceof Path ) {
            def converted = deepConvertToString(value)
            JsonOutput.prettyPrint(JsonOutput.toJson(converted))
        }

        else
            InvokerHelper.inspect(value)
    }

}
