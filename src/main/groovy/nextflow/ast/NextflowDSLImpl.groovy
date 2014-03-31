/*
 * Copyright (c) 2013-2014, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2014, Paolo Di Tommaso and the respective authors.
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
import nextflow.script.TokenEnvCall
import nextflow.script.TokenFileCall
import nextflow.script.TokenGString
import nextflow.script.TaskBody
import nextflow.script.TokenStdinCall
import nextflow.script.TokenStdoutCall
import nextflow.script.TokenValCall
import nextflow.script.TokenValRef
import nextflow.script.TokenVar
import org.codehaus.groovy.ast.ASTNode
import org.codehaus.groovy.ast.ClassCodeVisitorSupport
import org.codehaus.groovy.ast.ClassNode
import org.codehaus.groovy.ast.Parameter
import org.codehaus.groovy.ast.VariableScope
import org.codehaus.groovy.ast.expr.ArgumentListExpression
import org.codehaus.groovy.ast.expr.ClosureExpression
import org.codehaus.groovy.ast.expr.ConstantExpression
import org.codehaus.groovy.ast.expr.ConstructorCallExpression
import org.codehaus.groovy.ast.expr.Expression
import org.codehaus.groovy.ast.expr.GStringExpression
import org.codehaus.groovy.ast.expr.ListExpression
import org.codehaus.groovy.ast.expr.MapExpression
import org.codehaus.groovy.ast.expr.MethodCallExpression
import org.codehaus.groovy.ast.expr.PropertyExpression
import org.codehaus.groovy.ast.expr.TupleExpression
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
@CompileStatic
@GroovyASTTransformation(phase = CompilePhase.CONVERSION)
public class NextflowDSLImpl implements ASTTransformation {

    def String currentTaskName

    def String currentLabel

    @Override
    void visit(ASTNode[] astNodes, SourceUnit unit) {

        createVisitor(unit).visitClass((ClassNode)astNodes[1])

    }

    /*
     * create the code visitor
     */
    protected ClassCodeVisitorSupport createVisitor( SourceUnit unit ) {

        new ClassCodeVisitorSupport() {

            protected SourceUnit getSourceUnit() { unit }

            void visitMethodCallExpression(MethodCallExpression methodCall) {
                // pre-condition to be verified to apply the transformation
                Boolean preCondition = methodCall.objectExpression?.getText() == 'this'

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
    protected void convertProcessBlock( MethodCallExpression methodCall, SourceUnit unit ) {
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
            def source = new StringBuilder()
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
                            fixStdinStdout( stm )
                            convertInputMethod( stm.getExpression() )
                        }
                        break

                    case 'output':
                        if( stm instanceof ExpressionStatement ) {
                            fixStdinStdout( stm )
                            convertOutputMethod( stm.getExpression() )
                        }
                        break

                    case 'share':
                        if( stm instanceof ExpressionStatement ) {
                            convertShareMethod( stm.getExpression() )
                        }
                        break

                    case 'exec':
                        iterator.remove()
                        execStatements << stm
                        readSource(stm,source,unit)
                        break

                    case 'script':
                        iterator.remove()
                        execStatements << stm
                        readSource(stm,source,unit)
                        break

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
                // set the 'script' flag parameter
                def flag = currentLabel == 'script'
                def wrap = makeScriptWrapper(execClosure, source, flag, unit)
                block.addStatement( new ExpressionStatement(wrap)  )
                done = true

            }

            /*
             * when the last statement is a string script, the 'script:' label can be omitted
             */
            else if( len ) {
                def stm = block.getStatements().get(len-1)
                readSource(stm,source,unit)

                if ( stm instanceof ReturnStatement  ){
                    List result = wrapExpressionWithClosure(block, stm.getExpression(), len, source, unit)
                    done = result[0]
                    line = result[1] as int
                    coln = result[2] as int
                }

                else if ( stm instanceof ExpressionStatement )  {
                    List result = wrapExpressionWithClosure(block, stm.getExpression(), len, source, unit)
                    done = result[0]
                    line = result[1] as int
                    coln = result[2] as int
                }
                // set the 'script' flag
                currentLabel = 'script'
            }

            if (!done) {
                log.trace "Invalid 'process' definition -- Process must terminate with string expression"
                unit.addError( new SyntaxException("Not a valid process definition -- Make sure process ends with the script to be executed wrapped by quote characters", line,coln))
            }
        }
    }

    /**
     * Wrap the user provided piece of code, either a script or a closure with a {@code TaskBody} object
     *
     * @param closure
     * @param source
     * @param scriptOrNative
     * @param unit
     * @return a {@code TaskBody} object
     */
    private Expression makeScriptWrapper( ClosureExpression closure, CharSequence source, boolean scriptOrNative, SourceUnit unit ) {

        final newArgs = []
        newArgs << (closure)
        newArgs << ( new ConstantExpression(source.toString()) )
        newArgs << ( scriptOrNative ? ConstantExpression.PRIM_TRUE : ConstantExpression.PRIM_FALSE )

        final variables = fetchVariables(closure,unit)
        for( TokenValRef var: variables ) {
            def pName = new ConstantExpression(var.name)
            def pLine = new ConstantExpression(var.lineNum)
            def pCol = new ConstantExpression(var.colNum)
            newArgs << newObj( TokenValRef, pName, pLine, pCol )
        }

        newObj(TaskBody, newArgs as Object[] )
    }

    /**
     * Read the user provided script source string
     *
     * @param statement
     * @param buffer
     * @param unit
     */
    private void readSource( Statement statement, StringBuilder buffer, SourceUnit unit ) {

        def line = statement.getLineNumber()
        def last = statement.getLastLineNumber()
        for( int i=line; i<=last; i++ ) {
            buffer.append( unit.source.getLine(i, null) ) .append('\n')
        }

    }

    protected void fixStdinStdout( ExpressionStatement stm ) {

        if( stm.expression instanceof PropertyExpression ) {
            def expr = (PropertyExpression)stm.expression
            def obj = expr.objectExpression
            def prop = expr.property as ConstantExpression
            def target = new VariableExpression(prop.text)

            if( obj instanceof MethodCallExpression && 'stdout' == obj.getMethodAsString() ) {
                def stdout = new MethodCallExpression( new VariableExpression('this'), 'stdout', new ArgumentListExpression()  )
                def into = new MethodCallExpression(stdout, 'into', new ArgumentListExpression(target))
                // remove replace the old one with the new one
                stm.setExpression( into )
            }
            else if( obj instanceof MethodCallExpression && 'stdin' == obj.getMethodAsString() ) {
                def stdin = new MethodCallExpression( new VariableExpression('this'), 'stdin', new ArgumentListExpression()  )
                def from = new MethodCallExpression(stdin, 'from', new ArgumentListExpression(target))
                // remove replace the old one with the new one
                stm.setExpression( from )
            }
        }
    }

    /*
     * handle *input* parameters
     */
    protected void convertInputMethod( Expression expression ) {
        log.trace "convert > input expression: $expression"

        if( expression instanceof MethodCallExpression ) {

            def methodCall = expression as MethodCallExpression
            def methodName = methodCall.getMethodAsString()
            def nested = methodCall.objectExpression instanceof MethodCallExpression
            log.trace "convert > input method: $methodName"

            if( methodName in ['val','env','file','each', 'set','stdin'] ) {
                //this methods require a special prefix
                if( !nested )
                    methodCall.setMethod( new ConstantExpression('_in_' + methodName) )

                fixMethodCall(methodCall)
            }

            /*
             * Handles a GString a file name, like this:
             *
             *      input:
             *        file x name "$var_name" from q
             *
             */
            else if( methodName == 'name' && isWithinMethod(expression, 'file') ) {
                withinFileMethod = true
                varToConst(methodCall.getArguments())
                withinFileMethod = false
            }

            // invoke on the next method call
            if( expression.objectExpression instanceof MethodCallExpression ) {
                convertInputMethod(methodCall.objectExpression)
            }
        }

        else if( expression instanceof PropertyExpression ) {
            // invoke on the next method call
            if( expression.objectExpression instanceof MethodCallExpression ) {
                convertInputMethod(expression.objectExpression)
            }
        }

    }

    protected boolean isWithinMethod(MethodCallExpression method, String name) {
        if( method.objectExpression instanceof MethodCallExpression ) {
            return isWithinMethod(method.objectExpression as MethodCallExpression, name)
        }

        return method.getMethodAsString() == name
    }


    /*
     * handle *shared* parameters
     */

    protected void convertShareMethod( Expression expression ) {
        log.trace "convert > shared expression: $expression"

        if( expression instanceof MethodCallExpression ) {
            def methodCall = expression as MethodCallExpression
            def methodName = methodCall.getMethodAsString()
            def nested = methodCall.objectExpression instanceof MethodCallExpression
            log.trace "convert > shared method: $methodName"

            if( methodName in ['from','file','val','into','mode'] ) {
                if( !nested )
                    methodCall.setMethod( new ConstantExpression( '_share_' + methodName ) )
                fixMethodCall(methodCall)
            }

            if( methodCall.objectExpression instanceof MethodCallExpression ) {
                convertShareMethod(methodCall.objectExpression)
            }
        }
    }


    protected void convertOutputMethod( Expression expression ) {
        log.trace "convert > output expression: $expression"

        if( !(expression instanceof MethodCallExpression) ) {
            return
        }

        def methodCall = expression as MethodCallExpression
        def methodName = methodCall.getMethodAsString()
        def nested = methodCall.objectExpression instanceof MethodCallExpression
        log.trace "convert > output method: $methodName"

        if( methodName in ['val','file','set','stdout'] && !nested ) {
            // prefix the method name with the string '_out_'
            methodCall.setMethod( new ConstantExpression('_out_' + methodName) )
            fixMethodCall(methodCall)
        }

        else if( methodName in ['into','mode'] ) {
            fixMethodCall(methodCall)
        }

        // continue to traverse
        if( methodCall.objectExpression instanceof MethodCallExpression ) {
            convertOutputMethod(methodCall.objectExpression)
        }

    }

    private boolean withinSetMethod

    private boolean withinFileMethod

    /**
     * This method converts the a method call argument from a Variable to a Constant value
     * so that it is possible to reference variable that not yet exist
     *
     * @param methodCall The method object for which it is required to change args definition
     * @param flagVariable Whenever append a flag specified if the variable replacement has been applied
     * @param index The index of the argument to modify
     * @return
     */
    protected void fixMethodCall( MethodCallExpression methodCall ) {
        final name = methodCall.methodAsString

        withinSetMethod = name == '_in_set' || name == '_out_set'
        withinFileMethod = name == '_in_file' || name == '_out_file'

        try {
            varToConst(methodCall.getArguments())

        } finally {
            withinSetMethod = false
            withinFileMethod = false
        }

    }


    protected Expression varToStr( Expression expr ) {
        if( expr instanceof VariableExpression ) {
            def name = ((VariableExpression) expr).getName()
            return new ConstantExpression(name)
        }

        if( expr instanceof TupleExpression )  {
            def list = expr.getExpressions()
            list.eachWithIndex { Expression item, int i ->
                list[i] = varToStr(item)
            }
            return expr
        }

        return expr
    }

    protected Expression varToConst( Expression expr ) {

        if( expr instanceof VariableExpression ) {
            // when it is a variable expression, replace it with a constant representing
            // the variable name
            def name = ((VariableExpression) expr).getName()

            /*
             * the 'stdin' is used as placeholder for the standard input in the set definition. For example:
             *
             * input:
             *    set( stdin, .. ) from q
             */
            if( name == 'stdin' && withinSetMethod )
                return newObj( TokenStdinCall )

            /*
             * input:
             *    set( stdout, .. )
             */
            else if ( name == 'stdout' && withinSetMethod )
                return newObj( TokenStdoutCall )

            else
                return newObj( TokenVar, new ConstantExpression(name) )
        }

        /*
         * Handles GStrings in declaration file this
         * output:
         *   file "$name" into xx
         *   set( "$name" ) into xx
         *   set( file("$name") ) into xx
         */
        if( expr instanceof GStringExpression && ( withinFileMethod || withinSetMethod ) ) {

            def strings = (List<Expression>)expr.getStrings()
            def varNames = (List<Expression>)expr.getValues().collect { VariableExpression it -> new ConstantExpression(it.name) }

            return newObj( TokenGString, new ConstantExpression(expr.text), new ListExpression(strings), new ListExpression(varNames) )
        }

        /*
         * replace 'file' method call in the set definition, for example:
         *
         * input:
         *   set( file(fasta:'*.fa'), .. ) from q
         */
        if( expr instanceof MethodCallExpression && expr.methodAsString == 'file' && withinSetMethod ) {
            def args = (TupleExpression) varToConst(expr.arguments)
            return newObj( TokenFileCall, args )
        }

        /*
         * input:
         *  set( env(VAR_NAME) ) from q
         */
        if( expr instanceof MethodCallExpression && expr.methodAsString == 'env' && withinSetMethod ) {
            def args = (TupleExpression) varToStr(expr.arguments)
            return newObj( TokenEnvCall, args )
        }

        /*
         * input:
         *   set( val(x) ) from q
         */
        if( expr instanceof MethodCallExpression && expr.methodAsString == 'val' && withinSetMethod ) {
            def args = (TupleExpression) varToStr(expr.arguments)
            return newObj( TokenValCall, args )
        }

        // -- TupleExpression or ArgumentListExpression
        if( expr instanceof TupleExpression )  {
            def list = expr.getExpressions()
            list.eachWithIndex { Expression item, int i ->
                list[i] = varToConst(item)
            }
            return expr
        }

        return expr
    }

    /**
     * Creates a new {@code ConstructorCallExpression} for the specified class and arguments
     *
     * @param clazz The {@code Class} for which the create a constructor call expression
     * @param args The arguments to be passed to the constructor
     * @return The instance for the constructor call
     */
    protected Expression newObj( Class clazz, TupleExpression args ) {
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
    protected Expression newObj( Class clazz, Object... params ) {
        def type = new ClassNode(clazz)
        def args = new ArgumentListExpression( params as List<Expression>)
        return new ConstructorCallExpression(type,args)
    }

    /**
     * Wrap a generic expression with in a closure expression
     *
     * @param block The block to which the resulting closure has to be appended
     * @param expr The expression to the wrapped in a closure
     * @param len
     * @return A tuple in which:
     *      <li>1st item: {@code true} if successful or {@code false} otherwise
     *      <li>2nd item: on error condition the line containing the error in the source script, zero otherwise
     *      <li>3nd item: on error condition the column containing the error in the source script, zero otherwise
     *
     */
    protected List wrapExpressionWithClosure( BlockStatement block, Expression expr, int len, CharSequence source, SourceUnit unit ) {
        if( expr instanceof GStringExpression || expr instanceof ConstantExpression ) {
            // remove the last expression
            block.statements.remove(len-1)

            // and replace it by a wrapping closure
            def closureExp = new ClosureExpression( Parameter.EMPTY_ARRAY, new ExpressionStatement(expr) )
            closureExp.variableScope = new VariableScope(block.variableScope)

            // append to the list of statement
            //def wrap = newObj(TaskBody, closureExp, new ConstantExpression(source.toString()), ConstantExpression.TRUE)
            def wrap = makeScriptWrapper(closureExp, source, true, unit )
            block.statements.add( new ExpressionStatement(wrap) )

            return [true,0,0]
        }
        else if( expr instanceof ClosureExpression ) {
            // do not touch it
            return [true,0,0]
        }
        else {
            log.trace "Invalid process result expression: ${expr} -- Only constant or string expression can be used"
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
    protected void convertProcessDef( MethodCallExpression methodCall, SourceUnit unit ) {
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

    /**
     * Fetch all the variable references in a closure expression
     *
     * @param closure
     * @param unit
     * @return
     */
    protected Set<TokenValRef> fetchVariables( ClosureExpression closure, SourceUnit unit ) {
        def visitor = new VariableVisitor(unit)
        visitor.visitClosureExpression(closure)
        return visitor.allVariables
    }

    /**
     * Visit a closure and collect all referenced variable names
     */
    static class VariableVisitor extends ClassCodeVisitorSupport {

        final Map<String,TokenValRef> fAllVariables = [:]

        final SourceUnit sourceUnit

        VariableVisitor( SourceUnit unit ) {
            this.sourceUnit = unit
        }

        public void visitVariableExpression(VariableExpression var) {
            if( var.name == 'this' )
                return

            def token = new TokenValRef(var.name, var.lineNumber, var.columnNumber)
            NextflowDSLImpl.log.trace token.toString()

            if( !fAllVariables.containsKey(var.name))
                fAllVariables[var.name] = token
        }

        @Override
        protected SourceUnit getSourceUnit() {
            return sourceUnit
        }

        Set<TokenValRef> getAllVariables() {
            new HashSet<TokenValRef>(fAllVariables.values())
        }
    }


}
