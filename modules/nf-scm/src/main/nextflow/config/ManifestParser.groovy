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

package nextflow.config

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.config.ast.ConfigAssignNode
import nextflow.config.ast.ConfigBlockNode
import nextflow.config.ast.ConfigNode
import nextflow.config.ast.ConfigStatement
import nextflow.config.control.ConfigParser
import org.codehaus.groovy.ast.expr.ConstantExpression
import org.codehaus.groovy.ast.expr.Expression
import org.codehaus.groovy.ast.expr.ListExpression
import org.codehaus.groovy.ast.expr.MapExpression

/**
 * Reads the {@code manifest} scope of a Nextflow config file from its syntax tree,
 * without evaluating it, so repository code is never run to read a manifest.
 * Only literal values are read: strings, numbers, booleans, lists and maps.
 */
@Slf4j
@CompileStatic
class ManifestParser {

    private static final Object NOT_LITERAL = new Object()

    static Manifest parse(String text) {
        final source = new ConfigParser().parse('nextflow.config', text)
        if( source.getErrorCollector().hasErrors() )
            throw new IllegalArgumentException("Invalid config syntax")
        final result = new LinkedHashMap<String,Object>()
        collect((source.getAST() as ConfigNode).getConfigStatements(), [], result)
        return new Manifest(result)
    }

    private static void collect(List<ConfigStatement> statements, List<String> scope, Map<String,Object> result) {
        for( final stmt : statements ) {
            if( stmt instanceof ConfigBlockNode && stmt.kind == null ) {
                collect(stmt.statements, scope + stmt.name, result)
            }
            else if( stmt instanceof ConfigAssignNode ) {
                final names = scope + stmt.names
                if( names.size() != 2 || names[0] != 'manifest' )
                    continue
                final value = literal(stmt.value)
                if( value.is(NOT_LITERAL) )
                    log.debug "Ignoring non-literal manifest value: ${names[1]}"
                else
                    result.put(names[1], value)
            }
        }
    }

    private static Object literal(Expression expr) {
        if( expr instanceof ConstantExpression )
            return expr.getValue()
        if( expr instanceof ListExpression ) {
            final result = []
            for( final item : expr.getExpressions() ) {
                final value = literal(item)
                if( value.is(NOT_LITERAL) )
                    return NOT_LITERAL
                result.add(value)
            }
            return result
        }
        if( expr instanceof MapExpression ) {
            final result = new LinkedHashMap<Object,Object>()
            for( final entry : expr.getMapEntryExpressions() ) {
                final key = literal(entry.getKeyExpression())
                final value = literal(entry.getValueExpression())
                if( key.is(NOT_LITERAL) || value.is(NOT_LITERAL) )
                    return NOT_LITERAL
                result.put(key, value)
            }
            return result
        }
        return NOT_LITERAL
    }
}
