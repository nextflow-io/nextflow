/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import org.codehaus.groovy.ast.ASTNode
import org.codehaus.groovy.ast.ClassCodeExpressionTransformer
import org.codehaus.groovy.ast.ClassNode
import org.codehaus.groovy.ast.expr.BinaryExpression
import org.codehaus.groovy.ast.expr.ClosureExpression
import org.codehaus.groovy.ast.expr.Expression
import org.codehaus.groovy.ast.expr.MethodCallExpression
import org.codehaus.groovy.ast.expr.NotExpression
import org.codehaus.groovy.ast.tools.GeneralUtils
import org.codehaus.groovy.control.CompilePhase
import org.codehaus.groovy.control.SourceUnit
import org.codehaus.groovy.transform.ASTTransformation
import org.codehaus.groovy.transform.GroovyASTTransformation
/**
 * Implements Nextflow Xform logic
 * See http://groovy-lang.org/metaprogramming.html#_classcodeexpressiontransformer
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
@GroovyASTTransformation(phase = CompilePhase.CONVERSION)
class NextflowXformImpl implements ASTTransformation {

    SourceUnit unit

    @Override
    void visit(ASTNode[] nodes, SourceUnit source) {
        this.unit = unit
        createVisitor().visitClass((ClassNode)nodes[1])
    }

    protected ClassCodeExpressionTransformer createVisitor() {

        new ClassCodeExpressionTransformer() {

            protected SourceUnit getSourceUnit() { unit }

            @Override
            Expression transform(Expression expr) {
                if (expr == null)
                    return null

                def newExpr = transformBinaryExpression(expr)
                if( newExpr ) {
                    return newExpr
                }
                else if( expr instanceof ClosureExpression) {
                    visitClosureExpression(expr)
                }

                return super.transform(expr)
            }

            /**
             * This method replaces the `==` with the invocation of
             * {@link LangHelpers#compareEqual(java.lang.Object, java.lang.Object)}
             *
             * This is required to allow the comparisons of `Path` objects
             * which by default are not supported because it implements the Comparator interface
             *
             * See
             *  {@link LangHelpers#compareEqual(java.lang.Object, java.lang.Object)}
             *  https://stackoverflow.com/questions/28355773/in-groovy-why-does-the-behaviour-of-change-for-interfaces-extending-compar#comment45123447_28387391
             *
             */
            protected Expression transformBinaryExpression(Expression expr) {

                if( expr.class != BinaryExpression )
                    return null

                def binary = expr as BinaryExpression
                def left = binary.getLeftExpression()
                def right = binary.getRightExpression()

                if( '=='.equals(binary.operation.text) )
                    return call('compareEqual',left,right)

                if( '!='.equals(binary.operation.text) )
                    return new NotExpression(call('compareEqual',left,right))

                if( '<'.equals(binary.operation.text) )
                    return call('compareLessThan', left,right)

                if( '<='.equals(binary.operation.text) )
                    return call('compareLessThanEqual', left,right)

                if( '>'.equals(binary.operation.text) )
                    return call('compareGreaterThan', left,right)

                if( '>='.equals(binary.operation.text) )
                    return call('compareGreaterThanEqual', left,right)

                return null
            }


            private MethodCallExpression call(String method, Expression left, Expression right) {

                final a = transformBinaryExpression(left) ?: left
                final b = transformBinaryExpression(right) ?: right

                GeneralUtils.callX(
                        GeneralUtils.classX(LangHelpers),
                        method,
                        GeneralUtils.args(a,b))
            }

        }
    }

}
