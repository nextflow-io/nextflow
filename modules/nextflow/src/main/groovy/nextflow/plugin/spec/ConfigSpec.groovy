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
import nextflow.config.spec.SpecNode
import nextflow.script.types.Types
import org.codehaus.groovy.ast.ClassNode

/**
 * Generate specs for config scopes.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
class ConfigSpec {

    static Map<String,?> of(SpecNode node, String name) {
        return fromNode(node, name)
    }

    private static Map<String,?> fromNode(SpecNode node, String name) {
        if( node instanceof SpecNode.Option )
            return fromOption(node, name)
        if( node instanceof SpecNode.Placeholder )
            return fromPlaceholder(node, name)
        if( node instanceof SpecNode.Scope )
            return fromScope(node, name)
        throw new IllegalStateException()
    }

    private static Map<String,?> fromOption(SpecNode.Option node, String name) {
        final description = node.description().stripIndent(true).trim()
        final types = node.types().collect { t -> fromType(new ClassNode(t)) }

        return [
            type: 'ConfigOption',
            spec: [
                name: name,
                description: description,
                type: types.head(),
                additionalTypes: types.tail()
            ]
        ]
    }

    private static Map<String,?> fromPlaceholder(SpecNode.Placeholder node, String name) {
        final description = node.description().stripIndent(true).trim()
        final placeholderName = node.placeholderName()
        final scope = fromScope(node.scope())

        return [
            type: 'ConfigPlaceholderScope',
            spec: [
                name: name,
                description: description,
                placeholderName: placeholderName,
                scope: scope.spec
            ]
        ]
    }

    private static Map<String,?> fromScope(SpecNode.Scope node, String scopeName=null) {
        final description = node.description().stripIndent(true).trim()
        final children = node.children().collect { name, child ->
            fromNode(child, name)
        }

        return [
            type: 'ConfigScope',
            spec: [
                name: scopeName,
                description: description,
                children: children
            ]
        ]
    }

    private static Object fromType(ClassNode cn) {
        final name = Types.getName(cn.getTypeClass())
        if( !cn.isGenericsPlaceHolder() && cn.getGenericsTypes() != null ) {
            final typeArguments = cn.getGenericsTypes().collect { gt -> fromType(gt.getType()) }
            return [ name: name, typeArguments: typeArguments ]
        }
        else {
            return name
        }
    }
}
