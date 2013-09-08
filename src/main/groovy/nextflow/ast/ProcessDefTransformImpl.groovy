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

package nextflow.ast
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
class ProcessDefTransformImpl implements ASTTransformation {


    @Override
    void visit(ASTNode[] astNodes, SourceUnit unit) {

        createVisitor(unit).visitClass(astNodes[1])

    }

    /*
     * create the code visitor
     */
    def createVisitor( SourceUnit unit ) {

        new ClassCodeVisitorSupport() {

            def String currentTaskName

            protected SourceUnit getSourceUnit() { unit }

            void visitMethodCallExpression(MethodCallExpression methodCall) {

                // pre-condition to be verified to apply the transformation
                Boolean preCondition = methodCall.with {
                    (getMethod() instanceof ConstantExpression && objectExpression?.getText() == 'this')
                }

                /*
                 * intercept the *process* method in order to transform the script closure
                 */
                if( preCondition &&  methodCall.getMethodAsString() == 'process' ) {

                    currentTaskName = methodCall.getMethodAsString()
                    try {
                        convertProcessDef(methodCall,sourceUnit)
                        super.visitMethodCallExpression(methodCall)
                    }
                    finally {
                        currentTaskName = null
                    }

                }

                /*
                 * transform the 'input' method call invocation inside a 'task' method
                 */
                else if ( preCondition && methodCall.getMethodAsString() == 'input' && currentTaskName ) {

                    handleInputMethod(methodCall, sourceUnit)
                    super.visitMethodCallExpression(methodCall)
                }

                /*
                 * convert the *output* declarations
                 */
                else if ( preCondition && methodCall.getMethodAsString() == 'output' && currentTaskName  ) {
                    handleOutputMethod(methodCall, sourceUnit)
                    super.visitMethodCallExpression(methodCall)
                }

                /*
                 * convert the *stdout* declarations
                 */
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

        // need to convert the ArgumentListExpression to a 'TupleExpression'
        def args = methodCall.arguments as TupleExpression

        List<MapEntryExpression> entries = []
        args.getExpressions().each { Expression item ->

            /*
             * when the input declaration contains just a variable
             * it is converted to an 'val' input declaration  equivalent to
             *   val:'expression', from: expression
             */
            if ( item instanceof VariableExpression ) {
                entries << new MapEntryExpression( new ConstantExpression('val'), new ConstantExpression(item.getName()) )
                entries << new MapEntryExpression( new ConstantExpression('from'), item )
            }

            /*
             * convert input attributes from
             *    val: x
             * to
             *    val: 'x'
             *
             * so that users do not have to wrap values name by quotes
             */
            else if ( item instanceof MapExpression ) {
                item.mapEntryExpressions?.each { MapEntryExpression entry ->

                    def k = entry.keyExpression
                    def v = entry.valueExpression
                    def isConst = k instanceof ConstantExpression
                    if( isConst && k.value in ['val','each'] && v instanceof VariableExpression ) {
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

        log.trace "<< input converted args: ${entries}"
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

            /*
             * when the output declaration contains a constant value, it is
             * interpreted as the the file(s) name(s) to be returned by the task
             * the following map entry it is created
             *    file: 'name'
             */
            if( item instanceof ConstantExpression ) {
                entries << new MapEntryExpression( new ConstantExpression('file'), new ConstantExpression(item.getValue()) )
            }

            /*
             * when the output declaration contains just a variable value
             * it is interpreted as the the file(s) name(s) to be returned by the task
             * the following map entry it is created
             *    file: 'name'
             */
            else if ( item instanceof VariableExpression ) {
                entries << new MapEntryExpression( new ConstantExpression('file'), new ConstantExpression(item.getName()) )
            }

            else if ( item instanceof MapExpression ) {
                item.mapEntryExpressions?.each { MapEntryExpression entry ->
                    def k = entry.keyExpression
                    def v = entry.valueExpression

                    /*
                     * Channels in the output declarations have to referenced by name (string value).
                     * This fragment converts channel references from
                     *    into: variable
                     * to
                     *   into: 'name'
                     *
                     */
                    if( k instanceof ConstantExpression && k.value == 'into' && v instanceof VariableExpression ) {
                        entries << new MapEntryExpression( k, new ConstantExpression(v.name) )
                    }

                    /*
                     * Convert *val* reference like
                     *    val: variable
                     * to
                     *    val 'constant'
                     *
                     * so that the user do not have to wrap them by quotes
                     */
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

        log.trace "<< output converted args: ${entries}"
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

    /**
     * This method handle the process definition, so that it transform the user entered syntax
     *    process myName ( named: args, ..  ) { code .. }
     *
     * into
     *    process ( [named:args,..], String myName )  { }
     *
     * @param methodCall
     * @param unit
     */
    def void convertProcessDef( MethodCallExpression methodCall, SourceUnit unit ) {
        log.trace "Converts 'process' ${methodCall.arguments} "

        assert methodCall.arguments instanceof ArgumentListExpression
        def list = (methodCall.arguments as ArgumentListExpression).getExpressions()

        // extract the first argument which has to be a method-call expression
        // the name of this method represent the *process* name
        assert list.size() == 1
        assert list[0] instanceof MethodCallExpression
        def nested = list[0] as MethodCallExpression
        def name = nested.getMethodAsString()

        // the nested method arguments are the arguments to be passed
        // to the process definition, plus adding the process *name*
        // as an extra item in the arguments list
        def args = nested.getArguments() as ArgumentListExpression
        log.trace "Process name: $name with args: $args"

        // make sure to add the 'name' after the map item
        // (which represent the named parameter attributes)
        list = args.getExpressions()
        if( list.size()>0 && list[0] instanceof MapExpression ) {
            list.add(1, new ConstantExpression(name))
        }
        else {
            list.add(0, new ConstantExpression(name))
        }

        // set the new list as the new arguments
        methodCall.setArguments( args )

        // now continue as before !
        handleTaskMethod(methodCall, unit)
    }

}
