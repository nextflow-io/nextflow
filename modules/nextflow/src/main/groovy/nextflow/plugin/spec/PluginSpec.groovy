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
package nextflow.plugin.spec

import groovy.transform.CompileStatic
import nextflow.config.schema.ConfigScope
import nextflow.config.schema.SchemaNode
import nextflow.config.schema.ScopeName
import nextflow.plugin.extension.Factory
import nextflow.plugin.extension.Function
import nextflow.plugin.extension.Operator
import nextflow.plugin.extension.PluginExtensionPoint
import nextflow.script.dsl.Description

/**
 * Generate a plugin spec.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
class PluginSpec {

    private List<String> extensionPoints

    PluginSpec(List<String> extensionPoints) {
        this.extensionPoints = extensionPoints
    }

    Map build() {
        final classLoader = Thread.currentThread().getContextClassLoader()

        // extract schema for each plugin definition
        final definitions = []

        for( final className : extensionPoints ) {
            final clazz = classLoader.loadClass(className) as Class<? extends ConfigScope>

            if( ConfigScope.class.isAssignableFrom(clazz) ) {
                final scopeName = clazz.getAnnotation(ScopeName)?.value()
                final description = clazz.getAnnotation(Description)?.value()
                if( !scopeName )
                    continue
                final node = SchemaNode.Scope.of(clazz, description)

                definitions.add(ConfigSpec.of(node, scopeName))
            }

            if( PluginExtensionPoint.class.isAssignableFrom(clazz) ) {
                final methods = clazz.getDeclaredMethods()
                for( final method : methods ) {
                    if( method.getAnnotation(Factory) )
                        definitions.add(FunctionSpec.of(method, 'Factory'))
                    else if( method.getAnnotation(Function) )
                        definitions.add(FunctionSpec.of(method, 'Function'))
                    else if( method.getAnnotation(Operator) )
                        definitions.add(FunctionSpec.of(method, 'Operator'))
                }
            }
        }

        return [
            '$schema': 'https://raw.githubusercontent.com/nextflow-io/schemas/main/plugin/v1/schema.json',
            'definitions': definitions
        ]
    }
}
