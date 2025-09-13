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
package nextflow.plugin.index

import groovy.json.JsonOutput
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
 * Build an index of plugin definitions.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
class PluginIndex {

    private List<String> extensionPoints

    PluginIndex(List<String> extensionPoints) {
        this.extensionPoints = extensionPoints
    }

    List build() {
        final classLoader = Thread.currentThread().getContextClassLoader()

        // extract schema for each plugin definition
        final definitions = []

        for( final className : extensionPoints ) {
            final clazz = classLoader.loadClass(className)

            if( ConfigScope.class.isAssignableFrom(clazz) ) {
                final scopeName = clazz.getAnnotation(ScopeName)?.value()
                final description = clazz.getAnnotation(Description)?.value()
                if( !scopeName )
                    continue
                final node = SchemaNode.Scope.of(clazz, description)

                definitions.add(ConfigSchema.of(node, scopeName))
            }

            if( PluginExtensionPoint.class.isAssignableFrom(clazz) ) {
                final methods = clazz.getDeclaredMethods()
                for( final method : methods ) {
                    if( method.getAnnotation(Factory) )
                        definitions.add(FunctionSchema.of(method, 'Factory'))
                    else if( method.getAnnotation(Function) )
                        definitions.add(FunctionSchema.of(method, 'Function'))
                    else if( method.getAnnotation(Operator) )
                        definitions.add(FunctionSchema.of(method, 'Operator'))
                }
            }
        }

        return definitions
    }
}


@CompileStatic
class PluginIndexWriter {

    static void main(String[] args) {
        if (args.length < 2) {
            System.err.println("Usage: PluginIndexWriter <output-path> <class1> [class2] ...")
            System.exit(1)
        }

        final outputPath = args[0]
        final extensionPoints = args[1..-1]

        // build index
        final index = new PluginIndex(extensionPoints).build()

        // write index file
        final file = new File(outputPath)
        file.parentFile.mkdirs()
        file.text = JsonOutput.toJson(index)

        println "Saved plugin index to $file"
    }
}
