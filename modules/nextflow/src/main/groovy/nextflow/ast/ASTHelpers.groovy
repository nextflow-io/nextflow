/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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
import org.codehaus.groovy.ast.ASTNode
import org.codehaus.groovy.ast.ClassNode
import org.codehaus.groovy.ast.expr.ArgumentListExpression
import org.codehaus.groovy.ast.expr.BinaryExpression
import org.codehaus.groovy.ast.expr.CastExpression
import org.codehaus.groovy.ast.expr.ClosureExpression
import org.codehaus.groovy.ast.expr.ConstantExpression
import org.codehaus.groovy.ast.expr.ConstructorCallExpression
import org.codehaus.groovy.ast.expr.DeclarationExpression
import org.codehaus.groovy.ast.expr.Expression
import org.codehaus.groovy.ast.expr.MapExpression
import org.codehaus.groovy.ast.expr.MethodCallExpression
import org.codehaus.groovy.ast.expr.TupleExpression
import org.codehaus.groovy.ast.expr.VariableExpression
import org.codehaus.groovy.ast.stmt.BlockStatement
import org.codehaus.groovy.ast.stmt.ExpressionStatement
import org.codehaus.groovy.ast.stmt.ReturnStatement
import org.codehaus.groovy.ast.stmt.Statement
import org.codehaus.groovy.ast.tools.GeneralUtils
import org.codehaus.groovy.control.SourceUnit
import org.codehaus.groovy.syntax.SyntaxException
import org.codehaus.groovy.syntax.Types
/**
 * Helper methods to handle AST xforms
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class ASTHelpers {

    static void syntaxError(String message, ASTNode node, SourceUnit unit) {
        int line = node.lineNumber
        int coln = node.columnNumber
        final ex = new SyntaxException(message,line,coln)
        if( unit )
            unit.addError(ex)
        else
            throw ex
    }

    /**
     * Creates a new {@code ConstructorCallExpression} for the specified class and arguments
     *
     * @param clazz The {@code Class} for which the create a constructor call expression
     * @param args The arguments to be passed to the constructor
     * @return The instance for the constructor call
     */
    static Expression createX(Class clazz, TupleExpression args ) {
        def type = new ClassNode(clazz)
        return new ConstructorCallExpression(type,args)
    }

    /**
     * Creates a new {@code ConstructorCallExpression} for the specified class and arguments
     * specified using an open array. Te
     *
     * @param clazz The {@code Class} for which the create a constructor call expression
     * @param args The arguments to be passed to the constructor, they will be wrapped by in a {@code ArgumentListExpression}
     * @return The instance for the constructor call
     */
    static Expression createX(Class clazz, Expression... params ) {
        def type = new ClassNode(clazz)
        def args = new ArgumentListExpression(params as List<Expression>)
        return new ConstructorCallExpression(type,args)
    }

    static Expression createX(Class clazz, List<Expression> params ) {
        def type = new ClassNode(clazz)
        def args = new ArgumentListExpression(params as List<Expression>)
        return new ConstructorCallExpression(type,args)
    }

    static Expression declX(Expression left, Expression right) {
        new DeclarationExpression(left, GeneralUtils.ASSIGN, right)
    }

    static MethodCallExpression isMethodCallX(Expression expr) {
        return expr instanceof MethodCallExpression ? expr : null
    }

    static VariableExpression isVariableX(Expression expr) {
        return expr instanceof VariableExpression ? expr : null
    }

    static CastExpression isCastX(Expression expr) {
        return expr instanceof CastExpression ? expr : null
    }

    static VariableExpression isThisX(Expression expr) {
        isVariableX(expr)?.name == 'this' ? (VariableExpression)expr : null
    }

    static BinaryExpression isBinaryX(Expression expr) {
        expr instanceof BinaryExpression ? expr : null
    }

    static ArgumentListExpression isArgsX(Expression expr) {
        expr instanceof ArgumentListExpression ? expr : null
    }

    static TupleExpression isTupleX(Expression expr) {
        expr instanceof TupleExpression ? expr : null
    }

    static MapExpression isMapX(Expression exp) {
        exp instanceof MapExpression ? exp : null
    }

    static ConstantExpression isConstX(Expression exp) {
        exp instanceof ConstantExpression ? exp : null
    }

    static BinaryExpression isAssignX(Expression expr) {
        isBinaryX(expr)?.operation?.type == Types.ASSIGN ? (BinaryExpression)expr : null
    }

    static ClosureExpression isClosureX(Expression exp) {
        exp instanceof ClosureExpression ? exp : null
    }

    static BlockStatement isBlockStmt(Statement stm) {
        stm instanceof BlockStatement ? stm : null
    }

    static ExpressionStatement isStmtX(Statement stm) {
        stm instanceof ExpressionStatement ? stm : null
    }

    static ReturnStatement isReturnS(Statement stmt) {
        stmt instanceof ReturnStatement ? stmt : null
    }

}
