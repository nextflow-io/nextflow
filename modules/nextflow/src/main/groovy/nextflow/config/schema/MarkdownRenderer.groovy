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

import groovy.transform.TypeChecked
import nextflow.plugin.Plugins
import nextflow.script.dsl.Description

@TypeChecked
class MarkdownRenderer {

    String render() {
        final schema = getSchema()
        final entries = schema.entrySet().sort { entry -> entry.key }
        final result = new StringBuilder()
        entries.each { entry ->
            final scopeName = entry.key

            final anchor = scopeName == ''
                ? '(config-unscoped)='
                : "(config-$scopeName)="
            result.append("\n$anchor\n")

            final title = scopeName == ''
                ? 'Unscoped options'
                : "`$scopeName`"
            result.append("\n## $title\n")

            final scope = entry.value
            final description = scope.description()
            if( description )
                result.append("\n${fromDescription(description)}\n")
            result.append("\nThe following settings are available:\n")

            final options = scope.children().findAll { name, node -> node instanceof SchemaNode.Option }
            renderOptions(options, scopeName, result)

            final scopes = scope.children().findAll { name, node -> node instanceof SchemaNode.Scope }
            renderOptions(scopes, scopeName, result)
        }
        return result.toString()
    }

    private static Map<String,SchemaNode.Scope> getSchema() {
        final result = new HashMap<String,SchemaNode.Scope>()
        for( final scope : Plugins.getExtensions(ConfigScope) ) {
            final clazz = scope.getClass()
            final scopeName = clazz.getAnnotation(ScopeName)?.value()
            final description = clazz.getAnnotation(Description)?.value()
            if( scopeName == null )
                continue
            final node = SchemaNode.Scope.of(clazz, description)
            result.put(scopeName, node)
        }
        return result
    }

    private static String fromDescription(String description) {
        return description.stripIndent(true).trim()
    }

    private static void renderOptions(Map<String,SchemaNode> nodes, String scopeName, StringBuilder result) {
        final prefix = scopeName ? scopeName + '.' : ''
        final entries = nodes.entrySet().sort { entry -> entry.key }
        entries.each { entry ->
            final name = entry.key
            final node = entry.value
            if( node instanceof SchemaNode.Option )
                renderOption("${prefix}${name}", node, result)
            else if( node instanceof SchemaNode.Placeholder )
                renderOptions(node.scope().children(), "${prefix}${name}.${node.placeholderName()}", result)
            else if( node instanceof SchemaNode.Scope )
                renderOptions(node.children(), "${prefix}${name}", result)
            else
                throw new IllegalStateException()
        }
    }

    private static void renderOption(String name, SchemaNode.Option node, StringBuilder result) {
        final description = fromDescription(node.description())
        if( !description )
            return
        result.append("\n`${name}`")
        description.eachLine { line ->
            if( line )
                result.append("\n: ${line}")
        }
        result.append('\n')
    }

}
