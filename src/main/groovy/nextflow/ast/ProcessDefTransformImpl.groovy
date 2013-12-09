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
import org.codehaus.groovy.ast.expr.MapExpression
import org.codehaus.groovy.ast.expr.MethodCallExpression
import org.codehaus.groovy.ast.expr.VariableExpression
import org.codehaus.groovy.ast.stmt.BlockStatement
import org.codehaus.groovy.ast.stmt.ExpressionStatement
import org.codehaus.groovy.ast.stmt.ReturnStatement
import org.codehaus.groovy.ast.stmt.Statement
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

    def String currentTaskName

    def String currentLabel

    @Override
    void visit(ASTNode[] astNodes, SourceUnit unit) {

        createVisitor(unit).visitClass(astNodes[1])

    }

    /*
     * create the code visitor
     */
    def createVisitor( SourceUnit unit ) {

        new ClassCodeVisitorSupport() {


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

                    // clear block label
                    currentLabel = null
                    // keep track of 'process' method (it may be removed)
                    currentTaskName = methodCall.getMethodAsString()
                    try {
                        convertProcessDef(methodCall,sourceUnit)
                        super.visitMethodCallExpression(methodCall)
                    }
                    finally {
                        currentTaskName = null
                    }
                }

                // just apply the default behavior
                else {
                    super.visitMethodCallExpression(methodCall)
                }

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
    def void convertProcessBlock( MethodCallExpression methodCall, SourceUnit unit ) {
        log.trace "Apply task closure transformation to method call: $methodCall"

        def args = methodCall.arguments as ArgumentListExpression
        def lastArg = args.expressions.size()>0 ? args.getExpression(args.expressions.size()-1) : null
        def isClosure = lastArg instanceof ClosureExpression

        if( isClosure ) {
            // the block holding all the statements defined in the process (closure) definition
            def block = (lastArg as ClosureExpression).code as BlockStatement
            def len = block.statements.size()

            /*
             * iterate over the list of statements to:
             * - converts the method after the 'input:' label as input parameters
             * - converts the method after the 'output:' label as output parameters
             * - collect all the statement after the 'exec:' label
             */
            List<Statement> execStatements = []
            def iterator = block.getStatements().iterator()
            while( iterator.hasNext() ) {

                // get next statement
                Statement stm = iterator.next()

                // keep track of current block label
                currentLabel = stm.statementLabel ?: currentLabel

                switch(currentLabel) {
                    case 'input':
                        if( stm instanceof ExpressionStatement ) {
                            convertInputMethod( stm.getExpression() )
                        }
                        break

                    case 'output':
                        if( stm instanceof ExpressionStatement ) {
                            convertOutputMethod( stm.getExpression() )
                        }
                        break

                    case 'exec':
                        iterator.remove()
                        execStatements << stm
                }
            }

            /*
             * wrap all the statements after the 'exec:'  label by a new closure containing them (in a new block)
             */
            boolean done = false
            int line = methodCall.lineNumber
            int coln = methodCall.columnNumber
            if( execStatements ) {
                // create a new Closure
                def execBlock = new BlockStatement(execStatements, new VariableScope(block.variableScope))
                def execClosure = new ClosureExpression( Parameter.EMPTY_ARRAY, execBlock )

                // append the new block to the
                block.addStatement( new ExpressionStatement(execClosure) )
                done = true
            }

            /*
             * when the last statement is a string script, the 'exec:' label can be omitted
             */
            else if( len ) {
                def stm = block.getStatements().get(len-1)

                if ( stm instanceof ReturnStatement  ){
                    (done,line,coln) = wrapExpressionWithClosure(block, stm.expression, len)
                }

                else if ( stm instanceof ExpressionStatement )  {
                    (done,line,coln) = wrapExpressionWithClosure(block, stm.expression, len)
                }
            }


            if (!done) {
                log.trace "Not a valid task statement return type: ${stm.class.name} -- Task must terminate with string expression"
                unit.addError( new SyntaxException("Not a valid task definition -- Make sure task ends with the script to be executed wrapped by quote characters", line,coln))
            }
        }
    }

    def void convertInputMethod( Expression expression ) {
        log.trace "convert > input expression: $expression"

        if( !(expression instanceof MethodCallExpression) ) {
            return
        }

        def methodCall = expression as MethodCallExpression
        def methodName = methodCall.getMethodAsString()
        log.trace "convert > input method: $methodName"

        if( methodName in ['val','env','file','each'] ) {
            //this methods require a special prefix '__in_'
            methodCall.setMethod( new ConstantExpression('__in_' + methodName) )

            // the following methods require to replace a variable reference to a constant
            convertVarToConst(methodCall)

        }

        if( methodCall.objectExpression instanceof MethodCallExpression ) {
            convertInputMethod(methodCall.objectExpression)
        }

    }

    def void convertOutputMethod( Expression expression ) {
        log.trace "convert > output expression: $expression"

        if( !(expression instanceof MethodCallExpression) ) {
            return
        }

        def methodCall = expression as MethodCallExpression
        def methodName = methodCall.getMethodAsString()
        log.trace "convert > output method: $methodName"

        if( methodName in ['val','file'] ) {
            // prefix the method name with the string '__out_'
            methodCall.setMethod( new ConstantExpression('__out_' + methodName) )
        }

        if( methodName in ['val','file','using','stdout'] ) {
            convertVarToConst(methodCall)
        }

        // continue to traverse
        if( methodCall.objectExpression instanceof MethodCallExpression ) {
            convertOutputMethod(methodCall.objectExpression)
        }

    }

    /**
     * Converts a {@code VariableExpression} to a {@code ConstantExpression} having the same value as the variable name
     *
     * @param methodCall
     * @param index
     * @return
     */
    private List<Expression> convertVarToConst( MethodCallExpression methodCall, int index = 0 ) {

        def args = methodCall.getArguments() as ArgumentListExpression

        List<Expression> newArgs = []
        args.eachWithIndex { expr, i ->

            if( index == i && expr instanceof VariableExpression ) {
                newArgs << new ConstantExpression(expr.getName())
            }
            else {
                newArgs << expr
            }

        }
        methodCall.setArguments(new ArgumentListExpression(newArgs))

        return newArgs
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
        convertProcessBlock(methodCall, unit)
    }

}
