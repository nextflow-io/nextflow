/*
 * Copyright 2024-2025, Seqera Labs
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
package nextflow.config.schema

import groovy.json.JsonOutput
import groovy.transform.TypeChecked
import nextflow.plugin.Plugins
import nextflow.plugin.spec.ConfigSpec
import nextflow.script.dsl.Description

@TypeChecked
class JsonRenderer {

    String render() {
        final schema = getSchema()
        return JsonOutput.toJson(schema)
    }

    private static Map<String,?> getSchema() {
        final result = new HashMap<String,?>()
        for( final scope : Plugins.getExtensions(ConfigScope) ) {
            final clazz = scope.getClass()
            final scopeName = clazz.getAnnotation(ScopeName)?.value()
            final description = clazz.getAnnotation(Description)?.value()
            if( scopeName == '' ) {
                SchemaNode.Scope.of(clazz, '').children().each { name, node ->
                    result.put(name, ConfigSpec.of(node))
                }
                continue
            }
            if( !scopeName )
                continue
            final node = SchemaNode.Scope.of(clazz, description)
            result.put(scopeName, ConfigSpec.of(node, scopeName))
        }
        return result
    }

}
