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
 */

package nextflow.config

import java.nio.charset.StandardCharsets

import groovy.transform.CompileStatic
import org.codehaus.groovy.ast.ASTNode
import org.codehaus.groovy.ast.ClassCodeVisitorSupport
import org.codehaus.groovy.ast.ClassNode
import org.codehaus.groovy.ast.MethodNode
import org.codehaus.groovy.ast.expr.ConstantExpression
import org.codehaus.groovy.ast.stmt.BlockStatement
import org.codehaus.groovy.ast.stmt.ExpressionStatement
import org.codehaus.groovy.control.CompilePhase
import org.codehaus.groovy.control.SourceUnit
import org.codehaus.groovy.transform.ASTTransformation
import org.codehaus.groovy.transform.GroovyASTTransformation
/**
 * Implements Nextflow configuration file AST xform
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@GroovyASTTransformation(phase = CompilePhase.CONVERSION)
class ConfigTransformImpl implements ASTTransformation {

    @Override
    void visit(ASTNode[] astNodes, SourceUnit unit) {
        new ParserCodeVisitor(unit).visitClass((ClassNode)astNodes[1])
    }

    static class ParserCodeVisitor extends ClassCodeVisitorSupport {

        private SourceUnit unit

        ParserCodeVisitor(SourceUnit unit) {
            this.unit = unit
        }

        @Override
        protected SourceUnit getSourceUnit() { unit }

        @Override
        void visitMethod(MethodNode method) {
            if( method.name != 'run' )
                return

            // get embedded source code
            final block = (BlockStatement)method.code
            final stmt = (ExpressionStatement)block.statements.first()
            final strX = (ConstantExpression)stmt.expression
            final sourceBytes = ((String)strX.value).decodeBase64()
            final source = new String(sourceBytes, StandardCharsets.UTF_8)

            // build ast of embedded source using Nextflow parser
            final configUnit = ConfigAstBuilder.fromString(source)

            // insert ast
            method.code = new BlockStatement(
                configUnit.statements,
                block.variableScope
            )
        }

    }

}
