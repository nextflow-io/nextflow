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

package nextflow.extension

import java.nio.file.Path

import groovy.json.JsonOutput
import groovy.transform.CompileStatic
import nextflow.extension.FilesEx
import org.codehaus.groovy.runtime.InvokerHelper
import org.yaml.snakeyaml.DumperOptions
import org.yaml.snakeyaml.Yaml
/**
 * Utility functions for printing (dumping) values to standard output.
 *
 * @author : jorge <jorge.aguilera@seqera.io>
 */
@CompileStatic
class DumpHelper {

    /**
     * Normalize a value as follows:
     *
     * - file paths are serialized as URIs
     * - nulls are serialized using the given null value
     * - custom objects are serialized with toString()
     *
     * @param value
     * @param nullValue
     */
    static Object normalize(Object value, Object nullValue = null) {
        if( value instanceof Collection ) {
            return value.collect { it -> normalize(it) }
        }
        if( value instanceof Map ) {
            return value.collectEntries { k, v -> [k, normalize(v)] }
        }
        if( value instanceof Path ) {
            return FilesEx.toUriString(value)
        }
        if( value == null ) {
            return nullValue
        }
        return value.toString()
    }

    /**
     * Serialize a value as a pretty-printed string.
     *
     * @param value
     */
    static String prettyPrint(Object value) {
        if( value instanceof Collection || value instanceof Map || value instanceof Path )
            return prettyPrintJson(value)
        else
            return InvokerHelper.inspect(value)
    }

    /**
     * Serialize a value to a pretty-printed JSON string.
     *
     * @param value
     */
    static String prettyPrintJson(Object value) {
        return JsonOutput.prettyPrint(JsonOutput.toJson(normalize(value)))
    }

    /**
     * Serialize a value to a pretty-printed YAML string.
     *
     * Set `opts.style` to 'block' or 'flow' to optionally control
     * how arrays and objects are formatted.
     *
     * @param opts
     * @param value
     */
    static String prettyPrintYaml(Map opts = [:], Object value) {
        final dumperOptions = new DumperOptions()
        if( opts.style )
            dumperOptions.setDefaultFlowStyle(DumperOptions.FlowStyle.valueOf(opts.style.toString().toUpperCase()))
        dumperOptions.setSplitLines(false)
        return new Yaml(dumperOptions).dump(normalize(value))
    }

}
