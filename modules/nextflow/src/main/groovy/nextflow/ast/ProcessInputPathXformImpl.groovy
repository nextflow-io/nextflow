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

package nextflow.ast

import java.lang.reflect.Field
import java.lang.reflect.Modifier
import java.lang.reflect.ParameterizedType
import java.nio.file.Path

import static org.codehaus.groovy.ast.tools.GeneralUtils.*

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import org.codehaus.groovy.ast.ASTNode
import org.codehaus.groovy.ast.ClassCodeVisitorSupport
import org.codehaus.groovy.ast.ClassNode
import org.codehaus.groovy.ast.Parameter
import org.codehaus.groovy.ast.VariableScope
import org.codehaus.groovy.ast.expr.ArgumentListExpression
import org.codehaus.groovy.ast.expr.ClassExpression
import org.codehaus.groovy.ast.expr.ClosureExpression
import org.codehaus.groovy.ast.expr.ConstantExpression
import org.codehaus.groovy.ast.expr.Expression
import org.codehaus.groovy.ast.expr.MethodCallExpression
import org.codehaus.groovy.ast.expr.PropertyExpression
import org.codehaus.groovy.ast.expr.VariableExpression
import org.codehaus.groovy.ast.stmt.BlockStatement
import org.codehaus.groovy.ast.stmt.ExpressionStatement
import org.codehaus.groovy.ast.stmt.Statement
import org.codehaus.groovy.control.CompilePhase
import org.codehaus.groovy.control.SourceUnit
import org.codehaus.groovy.syntax.SyntaxException
import org.codehaus.groovy.transform.ASTTransformation
import org.codehaus.groovy.transform.GroovyASTTransformation

/**
 * Inject file staging directives for file inputs in
 * process definitions.
 *
 * Must be done during semantic analysis so that the
 * necessary type information is available.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
@GroovyASTTransformation(phase = CompilePhase.SEMANTIC_ANALYSIS)
class ProcessInputPathXformImpl implements ASTTransformation {

    @Override
    void visit(ASTNode[] astNodes, SourceUnit unit) {
        new DslCodeVisitor(unit).visitClass((ClassNode)astNodes[1])
    }

    @CompileStatic
    static class DslCodeVisitor extends ClassCodeVisitorSupport {

        private SourceUnit unit

        DslCodeVisitor(SourceUnit unit) {
            this.unit = unit
        }

        @Override
        protected SourceUnit getSourceUnit() { unit }

        @Override
        void visitMethodCallExpression(MethodCallExpression methodCall) {
            if( methodCall.objectExpression?.getText() != 'this' )
                return

            final methodName = methodCall.getMethodAsString()
            if( methodName != 'process' )
                return

            final args = methodCall.arguments as ArgumentListExpression
            final lastArg = args.expressions.size()>0 ? args.getExpression(args.expressions.size()-1) : null

            if( lastArg !instanceof ClosureExpression ) {
                syntaxError(lastArg, "Invalid process definition, possible syntax error")
                return
            }

            final closure = (ClosureExpression)lastArg
            final block = (BlockStatement)closure.code

            List<Statement> paramStatements = []
            String currentLabel = null

            for( final stmt : block.statements ) {
                currentLabel = stmt.statementLabel ?: currentLabel
                if( currentLabel != 'input' )
                    continue
                if( stmt !instanceof ExpressionStatement )
                    continue
                final stmtX = (ExpressionStatement)stmt
                if( stmtX.expression !instanceof MethodCallExpression )
                    continue
                final call = (MethodCallExpression)stmtX.expression
                if( call.methodAsString != '_typed_in_param' )
                    continue
                final paramArgs = (ArgumentListExpression)call.arguments
                assert paramArgs.size() == 2
                assert paramArgs[0] instanceof ConstantExpression
                assert paramArgs[1] instanceof ClassExpression
                final varName = ((ConstantExpression)paramArgs[0]).text
                final varType = ((ClassExpression)paramArgs[1]).type

                // infer staging directives via reflection
                final var = new VariableExpression(varName, varType)
                paramStatements.addAll( emitPathInputDecls(var, new TypeDef(var.type)) )
            }

            // prepend additional param statements
            paramStatements.addAll(block.statements)
            block.statements = paramStatements
        }

        protected List<Statement> emitPathInputDecls( Expression expr, TypeDef typeDef ) {
            List<Statement> result = []
            final type = typeDef.type

            if( isPathType(typeDef) ) {
                log.trace "inferring staging directive for path input: ${expr.text}"
                final block = new BlockStatement()
                block.addStatement( new ExpressionStatement(expr) )
                final closure = new ClosureExpression(Parameter.EMPTY_ARRAY, block)
                closure.variableScope = new VariableScope(block.variableScope)
                result << new ExpressionStatement( callThisX( 'stageAs', args(closure) ) )
            }
            else if( isRecordType(type) ) {
                for( final field : type.getDeclaredFields() ) {
                    if( /* !field.isAccessible() || */ Modifier.isStatic(field.getModifiers()) )
                        continue
                    log.trace "inspecting record type ${type.name}: field=${field.name}, type=${field.type.name}"
                    result.addAll( emitPathInputDecls(new PropertyExpression(expr, field.name), new TypeDef(field)) )
                }
            }

            return result
        }

        protected boolean isRecordType(Class type) {
            // NOTE: custom parser will be able to detect record types more elegantly
            log.trace "is ${type.name} a record type? ${type.package?.name == ''}"
            return type.package && type.package.name == ''
        }

        protected boolean isPathType(TypeDef typeDef) {
            final type = typeDef.type

            log.trace "is ${type.simpleName} a Path? ${type.name == 'java.nio.file.Path'}"
            if( Path.isAssignableFrom(type) )
                return true
            if( Collection.isAssignableFrom(type) && typeDef.genericTypes ) {
                final genericType = typeDef.genericTypes.first()
                log.trace "is ${type.simpleName}<${genericType.simpleName}> a Collection<Path>? ${genericType.name == 'java.nio.file.Path'}"
                return Path.isAssignableFrom(genericType)
            }
            return false
        }

        protected void syntaxError(ASTNode node, String message) {
            int line = node.lineNumber
            int coln = node.columnNumber
            unit.addError( new SyntaxException(message,line,coln))
        }

    }

    private static class TypeDef {
        Class type
        List<Class> genericTypes

        TypeDef(ClassNode classNode) {
            this.type = classNode.getPlainNodeReference().getTypeClass()
            if( classNode.getGenericsTypes() )
                this.genericTypes = classNode.getGenericsTypes().collect( el -> el.getType().getPlainNodeReference().getTypeClass() )
        }

        TypeDef(Field field) {
            this.type = field.type
            if( field.genericType instanceof ParameterizedType )
                this.genericTypes = field.genericType.getActualTypeArguments() as List<Class>
        }
    }

}
