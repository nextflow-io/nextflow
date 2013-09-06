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
                Boolean preCondition = methodCall.with {
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
                else if( preCondition &&  methodCall.getMethodAsString() == 'process' ) {

                    convertProcessDef(methodCall,sourceUnit)
                    super.visitMethodCallExpression(methodCall)
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

                else if ( preCondition && methodCall.getMethodAsString() == 'stdout' && currentTaskName  ) {
                    handleStdoutMethod(methodCall, sourceUnit)
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
     * Normalize the *input* declaration converting variables to fully qualified idiom
     * in the form {@code [ val: <literal>, from: <value> ] }
     *
     * @param methodCall
     * @param source
     */
    def void handleInputMethod( MethodCallExpression methodCall, SourceUnit source ) {
        log.trace "Method 'input' arguments: ${methodCall.arguments}"

        // if the arguments is already TupleExpression, there's nothing to do
        // since it is already wrapping a named parameters
        if ( methodCall.arguments.class == TupleExpression ) {
            log.trace "Transformation for argument of type: '${methodCall.arguments?.class?.name}' not required"
            return
        }

        // need to convert the ArgumentListExpression to a 'TupleExpression'
        def args = methodCall.arguments as ArgumentListExpression

        List<MapEntryExpression> entries = []
        args.getExpressions().each { Expression item ->
            if ( item instanceof VariableExpression ) {
                // when the input declaration contains just a variable
                // it is converted to an 'val' input declaration equivalent to
                // the following [ val:<literal>, channel:<value> ]
                entries << new MapEntryExpression( new ConstantExpression('val'), new ConstantExpression(item.getName()) )
                entries << new MapEntryExpression( new ConstantExpression('from'), item )
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

    /**
     * Normalize the {@code output} clause, interpreting constants and variables
     * as *output* files specification.
     *
     * @param methodCall
     * @param source
     */
    def void handleOutputMethod( MethodCallExpression methodCall, SourceUnit source ) {
        log.trace "Method 'output' arguments: ${methodCall.arguments}"

        // need to convert the ArgumentListExpression to a 'TupleExpression'
        def args = methodCall.arguments as TupleExpression

        List<MapEntryExpression> entries = []
        args.getExpressions().each { Expression item ->
            if( item instanceof ConstantExpression ) {
                // when the output declaration contains a constant value,
                // it is interpreted as the the file(s) name(s) to be returned by the task
                // the following map entry it is created {@code [ file:<const value> ] }
                entries << new MapEntryExpression( new ConstantExpression('file'), new ConstantExpression(item.getValue()) )
            }
            else if ( item instanceof VariableExpression ) {
                // when the output declaration contains just a variable value
                // it is interpreted as the the file(s) name(s) to be returned by the task
                // the following map entry it is created {@code [ file:<variable name> ] }
                entries << new MapEntryExpression( new ConstantExpression('file'), new ConstantExpression(item.getName()) )
            }

            else if ( item instanceof MapExpression ) {
                item.mapEntryExpressions?.each { MapEntryExpression entry ->
                    def k = entry.keyExpression
                    def v = entry.valueExpression

                    // converts {@code into: variable} to a {@code into: 'variable name' } clause
                    // in order to reference *channel* by name
                    if( k instanceof ConstantExpression && k.value == 'into' && v instanceof VariableExpression ) {
                        entries << new MapEntryExpression( k, new ConstantExpression(v.name) )
                    }

                    // converts val:x --> val:'x'
                    else if( k instanceof ConstantExpression && k.value == 'val' && v instanceof VariableExpression ) {
                        entries << new MapEntryExpression( k, new ConstantExpression(v.name) )
                    }

                    else {
                        entries << entry
                    }
                }
            }
            else {
                source.addError(new SyntaxException("Not a valid variable expression", item.lineNumber, item.columnNumber ))
                return
            }
        }

        TupleExpression tuple = new TupleExpression(new NamedArgumentListExpression(entries))
        methodCall.setArguments( tuple )

    }

    /**
     * Handle the 'stdout' task parameter definition. The 'stdout' it is supposed
     * to have an argument specified the 'channel' to which the task stdout has to be
     * forward.
     * <p>
     * Since we want to make possible to *create* a channel whenever it does not exist,
     * the channel has to be referenced by its name specified by a string value.
     * <p>
     * This method converts the clause {@code stdout channel} to {@code stdout 'channel'}
     *
     * @param methodCall
     * @param source
     */
    def void handleStdoutMethod( MethodCallExpression methodCall, SourceUnit source ) {
        log.trace "Method 'stdout' arguments: ${methodCall.arguments}"

        // need to convert the ArgumentListExpression to a 'TupleExpression'
        def args = methodCall.arguments as TupleExpression

        List<Expression> entries = []
        args.getExpressions().each { Expression item ->

            // converts variable expression
            // holding the reference to a channel
            // to a constant expression as he channel variable literal
            if ( item instanceof VariableExpression ) {
                entries << new ConstantExpression(item.getName())
            }
            else {
                source.addError(new SyntaxException("Not a valid stdout parameter value", item.lineNumber, item.columnNumber ))
                return
            }
        }

        TupleExpression tuple = new ArgumentListExpression(entries)
        methodCall.setArguments( tuple )

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
        log.trace "Apply task closure transformation to method call: $methodCall"

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


    def void convertProcessDef( MethodCallExpression methodCall, SourceUnit unit ) {
        log.debug "Converts 'process' ${methodCall.arguments} "

        assert methodCall.arguments instanceof ArgumentListExpression
        def list = (methodCall.arguments as ArgumentListExpression).getExpressions()

        assert list.size() == 1
        assert list[0] instanceof MethodCallExpression

        def nested = list[0] as MethodCallExpression
        def name = nested.getMethodAsString()

        def args = nested.getArguments() as ArgumentListExpression
        log.debug "Process name: $name with args: $args"

        args.getExpressions().add(0, new ConstantExpression(name))

        methodCall.setArguments( args )

        // now continue as before !
        handleTaskMethod(methodCall, unit)
    }

}
