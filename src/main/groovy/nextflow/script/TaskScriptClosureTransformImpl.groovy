/*
 * Copyright (c) 2012, the authors.
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

package nextflow.script
import groovy.util.logging.Slf4j
import org.codehaus.groovy.ast.ASTNode
import org.codehaus.groovy.ast.ClassCodeVisitorSupport
import org.codehaus.groovy.ast.Parameter
import org.codehaus.groovy.ast.VariableScope
import org.codehaus.groovy.ast.expr.ArgumentListExpression
import org.codehaus.groovy.ast.expr.ClosureExpression
import org.codehaus.groovy.ast.expr.ConstantExpression
import org.codehaus.groovy.ast.expr.Expression
import org.codehaus.groovy.ast.expr.GStringExpression
import org.codehaus.groovy.ast.expr.MapEntryExpression
import org.codehaus.groovy.ast.expr.MapExpression
import org.codehaus.groovy.ast.expr.MethodCallExpression
import org.codehaus.groovy.ast.expr.NamedArgumentListExpression
import org.codehaus.groovy.ast.expr.TupleExpression
import org.codehaus.groovy.ast.expr.VariableExpression
import org.codehaus.groovy.ast.stmt.BlockStatement
import org.codehaus.groovy.ast.stmt.ExpressionStatement
import org.codehaus.groovy.ast.stmt.ReturnStatement
import org.codehaus.groovy.control.CompilePhase
import org.codehaus.groovy.control.SourceUnit
import org.codehaus.groovy.syntax.SyntaxException
import org.codehaus.groovy.transform.ASTTransformation
import org.codehaus.groovy.transform.GroovyASTTransformation

/**
 * Implement some syntax sugars of Nextflow DSL scripting.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@Slf4j
@GroovyASTTransformation(phase = CompilePhase.CONVERSION)
class TaskScriptClosureTransformImpl implements ASTTransformation {


    @Override
    void visit(ASTNode[] nodes, SourceUnit unit) {

        createVisitor(unit).visitClass(nodes[1])

    }

    /**
     * Declares all method names to which is required to apply the closure transformation
     */
    static def final VALID_TASK_METHODS = ['task','merge']

    /*
     * create the code visitor
     */
    def createVisitor( SourceUnit unit ) {

        new ClassCodeVisitorSupport() {

            def String currentTaskName

            protected SourceUnit getSourceUnit() { unit }

            void visitMethodCallExpression(MethodCallExpression methodCall) {
                log.trace "Visiting methodCallExpr: ${methodCall}"

                // pre-condition to be verified to apply the transformation
                Boolean preCondition = methodCall.with{
                    (getMethod() instanceof ConstantExpression && objectExpression?.getText() == 'this')
                }

                // intercept the task method in order to transform the script closure
                if ( preCondition &&  methodCall.getMethodAsString() in VALID_TASK_METHODS ) {

                    currentTaskName = methodCall.getMethodAsString()
                    try {
                        handleTaskMethod(methodCall, sourceUnit)
                        super.visitMethodCallExpression(methodCall)
                    }
                    finally {
                        currentTaskName = null
                    }
                }

                // transform the 'input' method call invocation inside a 'task' method
                else if ( preCondition && methodCall.getMethodAsString() == 'input' && currentTaskName ) {

                    handleInputMethod(methodCall, sourceUnit)
                    super.visitMethodCallExpression(methodCall)
                }

                else if ( preCondition && methodCall.getMethodAsString() == 'output' && currentTaskName  ) {
                    handleOutputMethod(methodCall, sourceUnit)
                    super.visitMethodCallExpression(methodCall)
                }

                // just apply the default behavior
                else {
                    super.visitMethodCallExpression(methodCall)
                }

            }
        }
    }

    /**
     * This method transform method invocation for method {@code Processor#input}
     * from {@code input(x,y,z..) to {@code input( ['x':x, 'y':y, 'z':z ]}
     * <p>
     *     In other words from a list of values to a map for which each entry has the
     *     same name as the variable name itself
     *
     * @param methodCall
     * @param source
     */
    def void handleInputMethod( MethodCallExpression methodCall, SourceUnit source ) {
        log.trace "Method 'input' arguments: ${methodCall.arguments}"

        if ( !(methodCall.arguments instanceof ArgumentListExpression) ) {
            log.trace "Transformation for argument of type: '${methodCall.arguments?.class?.name}' not required"
            return
        }

        // need to convert the ArgumentListExpression to a 'TupleExpression'
        def args = methodCall.arguments as ArgumentListExpression

        List<MapEntryExpression> entries = []
        args.getExpressions().each { Expression item ->
            if ( item instanceof VariableExpression ) {
                   entries << new MapEntryExpression( new ConstantExpression(item.getName()), item )
            }
            else if ( item instanceof MapExpression ) {
                item.mapEntryExpressions?.each { MapEntryExpression entry -> entries << entry }
            }
            else {
                source.addError(new SyntaxException("Not a valid variable expression", item.lineNumber, item.columnNumber ))
                return
            }
        }

        TupleExpression tuple = new TupleExpression(new NamedArgumentListExpression(entries))
        methodCall.setArguments( tuple )

    }


    def void handleOutputMethod( MethodCallExpression methodCall, SourceUnit source ) {
        log.trace "Method 'output' arguments: ${methodCall.arguments}"

        if ( !(methodCall.arguments instanceof ArgumentListExpression) ) {
            log.trace "Transformation for argument of type: '${methodCall.arguments?.class?.name}' not required"
            return
        }

        // need to convert the ArgumentListExpression to a 'TupleExpression'
        def args = methodCall.arguments as ArgumentListExpression

        int c = args.getExpressions()?.size()
        for( int i=0; i<c; i++  ) {
            // replace any variable expression by a constant expression
            def item = args.expressions[i]
            if ( item instanceof VariableExpression ) {
                log.trace "Replacing Variable expression: $item"
                args.expressions[i] = new ConstantExpression(item.getName())
            }
            else {
                log.trace "Skipping expression item: $item"
            }
        }



    }

    /**
     * This transformation applies to tasks code block and translate
     * the script string to be executed wrapping it by a closure
     * <p>
     *     For example:
     * <pre>
     *     task {
     *         input x: someChannel
     *         output 'resultFile'
     *
     *         "do this; do that"
     *     }
     *
     * </pre>
     * becomes:
     *
     * <pre>
     *     task {
     *         input x: someChannel
     *         output 'resultFile'
     *
     *         { -> "do this; do that" }
     *     }
     *
     * </pre>
     *
     * @param methodCall
     * @param unit
     */
    def void handleTaskMethod( MethodCallExpression methodCall, SourceUnit unit ) {
        log.trace "Apply task closure trasformation to method call: $methodCall"

        def args = methodCall.arguments as ArgumentListExpression
        def lastArg = args.expressions.size()>0 ? args.getExpression(args.expressions.size()-1) : null
        def isClosure = lastArg instanceof ClosureExpression

        if( isClosure ) {
            def block = (lastArg as ClosureExpression).code as BlockStatement

            def len = block.statements.size()
            def stm = len>0 ? block.statements.get(len-1) : null

            boolean done = false
            int line = methodCall.lineNumber
            int coln = methodCall.columnNumber

            if ( stm instanceof ReturnStatement  ){
                (done,line,coln) = wrapExpressionWithClosure(block, stm.expression, len)
            }

            else if ( stm instanceof ExpressionStatement )  {
                (done,line,coln) = wrapExpressionWithClosure(block, stm.expression, len)
            }

            if (!done) {
                log.trace "Not a valid task statement return type: ${stm.class.name} -- Task must terminate with string expression"
                unit.addError( new SyntaxException("Not a valid task definition -- Make sure task ends with the script to be executed wrapped by quote characters", line,coln))
            }
        }
    }

    def List wrapExpressionWithClosure( BlockStatement block, Expression expr, int len ) {
        if( expr instanceof GStringExpression || expr instanceof ConstantExpression ) {
            // remove the last expression
            block.statements.remove(len-1)

            // and replace it by a wrapping closure
            def closureExp = new ClosureExpression( Parameter.EMPTY_ARRAY, new ExpressionStatement(expr) )
            closureExp.variableScope = new VariableScope()
            closureExp.variableScope.parent = block.variableScope

            // append to the list of statement
            block.statements.add( new ExpressionStatement(closureExp) )

            return [true,0,0]
        }
        else if( expr instanceof ClosureExpression ) {
            // do not touch it
            return [true,0,0]
        }
        else {
            log.trace "Invalid task result expression: ${expr} -- Only constant or string expression can be used"
        }

        return [false, expr.lineNumber, expr.columnNumber]
    }



}
