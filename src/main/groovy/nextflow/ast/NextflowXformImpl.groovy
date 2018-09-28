/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
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
                if( expr.class == BinaryExpression ) {
                    def newExpr = replaceBinaryExpression(expr as BinaryExpression)
                    if( newExpr )
                        return newExpr
                }
                else if( expr instanceof ClosureExpression) {
                    visitClosureExpression(expr)
                }

                super.transform(expr)
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
            protected Expression replaceBinaryExpression(BinaryExpression expr) {

                final binary = expr as BinaryExpression
                final left = binary.getLeftExpression()
                final right = binary.getRightExpression()

                if( '=='.equals(binary.operation.text) ) {
                    return call('compareEqual',left,right)
                }

                if( '!='.equals(binary.operation.text) ) {
                    return new NotExpression(call('compareEqual',left,right))
                }

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


            private static MethodCallExpression call(String method, Expression left, Expression right) {
                GeneralUtils.callX(
                        GeneralUtils.classX(LangHelpers),
                        method,
                        GeneralUtils.args(left,right))
            }

        }
    }

}
