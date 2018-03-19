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

package nextflow.config

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import org.codehaus.groovy.ast.ASTNode
import org.codehaus.groovy.ast.AnnotationNode
import org.codehaus.groovy.ast.ClassCodeVisitorSupport
import org.codehaus.groovy.ast.ClassNode
import org.codehaus.groovy.ast.expr.ArgumentListExpression
import org.codehaus.groovy.ast.expr.BinaryExpression
import org.codehaus.groovy.ast.expr.ClosureExpression
import org.codehaus.groovy.ast.expr.ConstantExpression
import org.codehaus.groovy.ast.expr.ConstructorCallExpression
import org.codehaus.groovy.ast.expr.Expression
import org.codehaus.groovy.ast.expr.MapEntryExpression
import org.codehaus.groovy.ast.expr.MapExpression
import org.codehaus.groovy.ast.expr.MethodCallExpression
import org.codehaus.groovy.ast.stmt.ExpressionStatement
import org.codehaus.groovy.control.CompilePhase
import org.codehaus.groovy.control.SourceUnit
import org.codehaus.groovy.transform.ASTTransformation
import org.codehaus.groovy.transform.GroovyASTTransformation
/**
 * Implements Nextflow configuration file AST xform
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
@GroovyASTTransformation(phase = CompilePhase.CONVERSION)
class ConfigTransformImpl implements ASTTransformation {

    private boolean renderClosureAsString

    @Override
    void visit(ASTNode[] astNodes, SourceUnit unit) {
        final annot = (AnnotationNode)astNodes[0]
        final clazz = (ClassNode)astNodes[1]
        // the following line is mostly an hack to pass a parameter to this xform instance
        this.renderClosureAsString = annot.getMember('renderClosureAsString') != null
        createVisitor(unit).visitClass(clazz)
    }

    protected ClassCodeVisitorSupport createVisitor(SourceUnit unit) {
        return (renderClosureAsString
                ? new RetainClosureSourceCodeVisitorSupport(unit: unit)
                : new DefaultConfigCodeVisitor(unit: unit) )
    }

    /**
     * Nextflow config file visitor support class. Apply default transformations to
     * the config object to implements config syntax sugars
     */
    @CompileStatic
    static class DefaultConfigCodeVisitor extends ClassCodeVisitorSupport {

        protected SourceUnit unit

        @Override
        protected SourceUnit getSourceUnit() { unit }

        @Override
        void visitExpressionStatement(ExpressionStatement stm) {
            if( stm.expression instanceof MethodCallExpression && stm.getStatementLabel() == 'withLabel' ) {
                replaceMethodName( stm.expression as MethodCallExpression, 'withLabel' )
            }
            else if( stm.expression instanceof MethodCallExpression && stm.getStatementLabel() == 'withName' ) {
                replaceMethodName( stm.expression as MethodCallExpression, 'withName' )
            }
            super.visitExpressionStatement(stm)
        }

        /**
         * Replace the name of the invoked method pre-pending with the specified string
         *
         * @param call A object representing a method call expression
         * @param prefix A string to prepend to the name of the invoked method
         */
        protected void replaceMethodName(MethodCallExpression call, String prefix) {
            call.setMethod( new ConstantExpression(prefix + ":" + call.method.text) )
        }

    }

    /**
     * This visitor is only used to render the closure source code
     * when is required to visualise the nextflow config content
     */
    @CompileStatic
    static class RetainClosureSourceCodeVisitorSupport extends DefaultConfigCodeVisitor {

        /**
         * Visit expression statements replacing a binary assignment such as:
         *
         *      foo = { closure }
         *
         * with an equivalent assignment replacing the closure with placeholder
         * class holding the closure source code
         *
         * @param stm
         *
         * @see ConfigClosurePlaceholder
         */
        @Override
        void visitExpressionStatement(ExpressionStatement stm) {
            if( stm.expression instanceof BinaryExpression ) {
                replaceClosureAssignment(stm, stm.expression as BinaryExpression)
            }
            super.visitExpressionStatement(stm)
        }

        protected void replaceClosureAssignment(ExpressionStatement stm, BinaryExpression expr ) {
            if( expr.operation.text == '=' && expr.rightExpression instanceof ClosureExpression ) {
                final value = closureToString(expr.rightExpression)
                final replace = new BinaryExpression(expr.leftExpression, expr.getOperation(), value)
                stm.setExpression(replace)
            }
        }

        /**
         * Visit Map expressions replacing closure values such as :
         *
         *      [key: { closure }]
         *
         * with an equivalent Map in which the closure is replaced with a placeholder
         * class holding the closure source code
         *
         * @param expr
         *
         * @see ConfigClosurePlaceholder
         */
        @Override
        void visitMapExpression(MapExpression expr) {
            for( int i=0; i<expr.mapEntryExpressions.size(); i++ ) {
                def entry = expr.mapEntryExpressions[i]
                if( entry.valueExpression instanceof ClosureExpression ) {
                    expr.mapEntryExpressions[i] = replaceMapEntryAssignment(entry)
                }
            }
            super.visitMapExpression(expr)
        }

        protected MapEntryExpression replaceMapEntryAssignment(MapEntryExpression expr) {
            if( expr.valueExpression instanceof ClosureExpression ) {
                final value = closureToString(expr.valueExpression)
                new MapEntryExpression(expr.keyExpression, value)
            }
            else
                return expr
        }

        protected Expression closureToString( Expression closure )  {
            def buffer = new StringBuilder()
            readSource(closure, buffer)
            def str = new ConstantExpression(buffer.toString())

            def type = new ClassNode(ConfigClosurePlaceholder)
            def args = new ArgumentListExpression( [str] as List<Expression> )
            return new ConstructorCallExpression(type,args)
        }

        /**
         * Read the user provided script source string
         *
         * @param expr
         * @param buffer
         * @param unit
         */
        protected void readSource(Expression expr, StringBuilder buffer) {
            final colBegin = Math.max(expr.getColumnNumber()-1, 0)
            final colEnd = Math.max(expr.getLastColumnNumber()-1, 0)
            final lineFirst = expr.getLineNumber()
            final lineLast = expr.getLastLineNumber()

            for( int i=lineFirst; i<=lineLast; i++ ) {
                def line = unit.source.getLine(i, null)
                if( i==lineFirst ) {
                    def str = i==lineLast ? line.substring(colBegin,colEnd) : line.substring(colBegin)
                    buffer.append(str)
                }
                else {
                    def str = i==lineLast ? line.substring(0, colEnd) : line
                    buffer.append('\n')
                    buffer.append(str)
                }
            }
        }

    }

}
