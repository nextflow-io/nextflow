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
import nextflow.script.dsl.Description
import nextflow.script.types.Types
import org.codehaus.groovy.ast.ClassNode

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
                    result.put(name, fromNode(node))
                }
                continue
            }
            if( !scopeName )
                continue
            final node = SchemaNode.Scope.of(clazz, description)
            result.put(scopeName, fromNode(node, scopeName))
        }
        return result
    }

    private static Map<String,?> fromNode(SchemaNode node, String name=null) {
        if( node instanceof SchemaNode.Option )
            return fromOption(node)
        if( node instanceof SchemaNode.Placeholder )
            return fromPlaceholder(node)
        if( node instanceof SchemaNode.Scope )
            return fromScope(node, name)
        throw new IllegalStateException()
    }

    private static Map<String,?> fromOption(SchemaNode.Option node) {
        final description = node.description().stripIndent(true).trim()
        final type = fromType(new ClassNode(node.type()))

        return [
            type: 'Option',
            spec: [
                description: description,
                type: type
            ]
        ]
    }

    private static Map<String,?> fromPlaceholder(SchemaNode.Placeholder node) {
        final description = node.description().stripIndent(true).trim()
        final placeholderName = node.placeholderName()
        final scope = fromScope(node.scope())

        return [
            type: 'Placeholder',
            spec: [
                description: description,
                placeholderName: placeholderName,
                scope: scope.spec
            ]
        ]
    }

    private static Map<String,?> fromScope(SchemaNode.Scope node, String scopeName=null) {
        final description = node.description().stripIndent(true).trim()
        final children = node.children().collectEntries { name, child ->
            Map.entry(name, fromNode(child, name))
        }

        return [
            type: 'Scope',
            spec: [
                description: withLink(scopeName, description),
                children: children
            ]
        ]
    }

    private static String withLink(String scopeName, String description) {
        return scopeName
            ? "$description\n\n[Read more](https://nextflow.io/docs/latest/reference/config.html#$scopeName)\n"
            : description
    }

    private static Object fromType(ClassNode cn) {
        final name = Types.getName(cn.getTypeClass())
        if( cn.isUsingGenerics() ) {
            final typeArguments = cn.getGenericsTypes().collect { gt -> fromType(gt.getType()) }
            return [ name: name, typeArguments: typeArguments ]
        }
        else {
            return name
        }
    }

}
