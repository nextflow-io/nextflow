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
import nextflow.script.BaseScript
import nextflow.script.TaskBody
import nextflow.script.TaskClosure
import nextflow.script.TokenEnvCall
import nextflow.script.TokenFileCall
import nextflow.script.TokenStdinCall
import nextflow.script.TokenStdoutCall
import nextflow.script.TokenValCall
import nextflow.script.TokenValRef
import nextflow.script.TokenVar
import org.codehaus.groovy.ast.ASTNode
import org.codehaus.groovy.ast.ClassCodeVisitorSupport
import org.codehaus.groovy.ast.ClassNode
import org.codehaus.groovy.ast.ConstructorNode
import org.codehaus.groovy.ast.Parameter
import org.codehaus.groovy.ast.VariableScope
import org.codehaus.groovy.ast.expr.ArgumentListExpression
import org.codehaus.groovy.ast.expr.BinaryExpression
import org.codehaus.groovy.ast.expr.ClosureExpression
import org.codehaus.groovy.ast.expr.ConstantExpression
import org.codehaus.groovy.ast.expr.ConstructorCallExpression
import org.codehaus.groovy.ast.expr.DeclarationExpression
import org.codehaus.groovy.ast.expr.Expression
import org.codehaus.groovy.ast.expr.GStringExpression
import org.codehaus.groovy.ast.expr.ListExpression
import org.codehaus.groovy.ast.expr.MapExpression
import org.codehaus.groovy.ast.expr.MethodCallExpression
import org.codehaus.groovy.ast.expr.PropertyExpression
import org.codehaus.groovy.ast.expr.TupleExpression
import org.codehaus.groovy.ast.expr.UnaryMinusExpression
import org.codehaus.groovy.ast.expr.VariableExpression
import org.codehaus.groovy.ast.stmt.BlockStatement
import org.codehaus.groovy.ast.stmt.ExpressionStatement
import org.codehaus.groovy.ast.stmt.ReturnStatement
import org.codehaus.groovy.ast.stmt.Statement
import org.codehaus.groovy.control.CompilePhase
import org.codehaus.groovy.control.SourceUnit
import org.codehaus.groovy.syntax.SyntaxException
import org.codehaus.groovy.syntax.Types
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

    protected String currentTaskName

    protected String currentLabel

    protected GStringToLazyVisitor makeGStringLazyVisitor

    protected Set<String> processNames = []

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

            /**
             * Creates a statement that invokes the {@link nextflow.script.BaseScript#init(java.util.List)} method
             * used to initialize the script with metadata collected during script parsing
             * @return The method invocation statement
             */
            protected Statement invokeBaseScriptInit() {
                final names = new ListExpression()
                processNames.each { names.addExpression(new ConstantExpression(it.toString())) }

                // the method list argument
                final args = new ArgumentListExpression()
                args.addExpression(names)

                final call = new MethodCallExpression(new VariableExpression('this'), 'init', args)
                final stm = new ExpressionStatement(call)
                return stm
            }

            /**
             * Add to constructor a method call to inject parsed metadata
             * @param node
             */
            protected void injectMetadata(ClassNode node) {
                for( ConstructorNode constructor : node.getDeclaredConstructors() ) {
                    def code = constructor.getCode()
                    if( code instanceof BlockStatement ) {
                        code.addStatement( invokeBaseScriptInit() )
                    }
                    else if( code instanceof ExpressionStatement ) {
                        def expr = code
                        def block = new BlockStatement()
                        block.addStatement(expr)
                        block.addStatement( invokeBaseScriptInit() )
                        constructor.setCode(block)
                    }
                    else
                        throw new IllegalStateException("Invalid constructor expression: $code")
                }
            }

            @Override
            protected void visitObjectInitializerStatements(ClassNode node) {
                if( node.getSuperClass().getName() == BaseScript.getName() ) {
                    injectMetadata(node)
                }
                super.visitObjectInitializerStatements(node)
            }

            @Override
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
     * Transform a DSL `process` definition into a proper method invocation
     *
     * @param methodCall
     * @param unit
     */
    protected void convertProcessBlock( MethodCallExpression methodCall, SourceUnit unit ) {
        log.trace "Apply task closure transformation to method call: $methodCall"

        final args = methodCall.arguments as ArgumentListExpression
        final lastArg = args.expressions.size()>0 ? args.getExpression(args.expressions.size()-1) : null
        final isClosure = lastArg instanceof ClosureExpression

        if( isClosure ) {
            // the block holding all the statements defined in the process (closure) definition
            final block = (lastArg as ClosureExpression).code as BlockStatement

            makeGStringLazyVisitor = new GStringToLazyVisitor(unit)

            /*
             * iterate over the list of statements to:
             * - converts the method after the 'input:' label as input parameters
             * - converts the method after the 'output:' label as output parameters
             * - collect all the statement after the 'exec:' label
             */
            def source = new StringBuilder()
            List<Statement> execStatements = []

            List<Statement> whenStatements = []
            def whenSource = new StringBuilder()

            def iterator = block.getStatements().iterator()
            while( iterator.hasNext() ) {

                // get next statement
                Statement stm = iterator.next()

                // keep track of current block label
                currentLabel = stm.statementLabel ?: currentLabel

                switch(currentLabel) {
                    case 'input':
                        if( stm instanceof ExpressionStatement ) {
                            fixLazyGString( stm )
                            fixStdinStdout( stm )
                            convertInputMethod( stm.getExpression() )
                        }
                        break

                    case 'output':
                        if( stm instanceof ExpressionStatement ) {
                            fixLazyGString( stm )
                            fixStdinStdout( stm )
                            convertOutputMethod( stm.getExpression() )
                        }
                        break

                    case 'exec':
                        iterator.remove()
                        execStatements << stm
                        readSource(stm,source,unit)
                        break

                    case 'script':
                    case 'shell':
                        iterator.remove()
                        execStatements << stm
                        readSource(stm,source,unit)
                        break

                    // capture the statements in a when guard and remove from the current block
                    case 'when':
                        if( iterator.hasNext() ) {
                            iterator.remove()
                            whenStatements << stm
                            readSource(stm,whenSource,unit)
                            break
                        }
                        // when entering in this branch means that this is the last statement,
                        // which is supposed to be the task command
                        // hence if no previous `when` statement has been processed, a syntax error is returned
                        else if( !whenStatements ) {
                            int line = methodCall.lineNumber
                            int coln = methodCall.columnNumber
                            unit.addError(new SyntaxException("Invalid process definition -- Empty `when` or missing `script` statement", line, coln))
                            return
                        }
                        else
                            break

                    default:
                        if(currentLabel) {
                            def line = stm.getLineNumber()
                            def coln = stm.getColumnNumber()
                            unit.addError(new SyntaxException("Invalid process definition -- Unknown keyword `$currentLabel`",line,coln))
                            return
                        }

                        fixLazyGString(stm)
                        fixDirectiveWithNegativeValue(stm)  // Fixes #180
                }
            }

            /*
             * add the `when` block if found
             */
            if( whenStatements ) {
                addWhenGuardCall(whenStatements, whenSource, block)
            }

            /*
             * wrap all the statements after the 'exec:'  label by a new closure containing them (in a new block)
             */
            final len = block.statements.size()
            boolean done = false
            if( execStatements ) {
                // create a new Closure
                def execBlock = new BlockStatement(execStatements, new VariableScope(block.variableScope))
                def execClosure = new ClosureExpression( Parameter.EMPTY_ARRAY, execBlock )

                // append the new block to the
                // set the 'script' flag parameter
                def wrap = makeScriptWrapper(execClosure, source, currentLabel, unit)
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
                    done = wrapExpressionWithClosure(block, stm.getExpression(), len, source, unit)
                }

                else if ( stm instanceof ExpressionStatement )  {
                    done = wrapExpressionWithClosure(block, stm.getExpression(), len, source, unit)
                }

                // set the 'script' flag
                currentLabel = 'script'
            }

            if (!done) {
                log.trace "Invalid 'process' definition -- Process must terminate with string expression"
                int line = methodCall.lineNumber
                int coln = methodCall.columnNumber
                unit.addError( new SyntaxException("Invalid process definition -- Make sure the process ends with a script wrapped by quote characters",line,coln))
            }
        }
    }

    /**
     * Converts a `when` block into a when method call expression. The when code is converted into a
     * closure expression and set a `when` directive in the process configuration properties.
     *
     * See {@link nextflow.processor.ProcessConfig#configProperties}
     * See {@link nextflow.processor.TaskConfig#getGuard(java.lang.String)}
     */
    protected void addWhenGuardCall( List<Statement> statements, StringBuilder source, BlockStatement parent ) {
        // wrap the code block into a closure expression
        def block = new BlockStatement(statements, new VariableScope(parent.variableScope))
        def closure = new ClosureExpression( Parameter.EMPTY_ARRAY, block )

        // the closure expression is wrapped itself into a TaskClosure object
        // in order to capture the closure source other than the closure code
        def newArgs = []
        newArgs << closure
        newArgs << new ConstantExpression(source.toString())
        def whenObj = newObj( TaskClosure, newArgs as Object[] )

        // creates a method call expression for the method `when`
        def method = new MethodCallExpression(VariableExpression.THIS_EXPRESSION, 'when', whenObj)
        parent.getStatements().add(0, new ExpressionStatement(method))

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
    private Expression makeScriptWrapper( ClosureExpression closure, CharSequence source, String section, SourceUnit unit ) {

        final newArgs = []
        newArgs << (closure)
        newArgs << ( new ConstantExpression(source.toString()) )
        newArgs << ( new ConstantExpression(section) )

        final variables = fetchVariables(closure,unit)
        for( TokenValRef var: variables ) {
            def pName = new ConstantExpression(var.name)
            def pLine = new ConstantExpression(var.lineNum)
            def pCol = new ConstantExpression(var.colNum)
            newArgs << newObj( TokenValRef, pName, pLine, pCol )
        }

        // invokes the TaskBody constructor
        newObj( TaskBody, newArgs as Object[] )
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

    @Deprecated
    protected void fixGuards( Statement stm, SourceUnit unit ) {

        if( stm instanceof ExpressionStatement && stm.getExpression() instanceof MethodCallExpression ) {

            def method = stm.getExpression() as MethodCallExpression
            if( method.arguments instanceof ArgumentListExpression ) {
                def args = method.arguments as ArgumentListExpression
                if( args.size() == 1 && args.getExpression(0) instanceof ClosureExpression ) {
                    // read the source code of the closure
                    def closure = args.getExpression(0) as ClosureExpression
                    def source = new StringBuilder()
                    readSource(closure.getCode(), source, unit)

                    // wrap the closure expression by a TaskClosure object invocation
                    def wrap = newObj( TaskClosure, closure, new ConstantExpression(source.toString()) )
                    // replace it with the original closure argument
                    args.expressions.set(0, wrap)
                }
            }
        }

    }

    protected void fixLazyGString( Statement stm ) {
        if( stm instanceof ExpressionStatement && stm.getExpression() instanceof MethodCallExpression ) {
            makeGStringLazyVisitor.visitExpressionStatement(stm)
        }
    }

    protected void fixDirectiveWithNegativeValue( Statement stm ) {
        if( stm instanceof ExpressionStatement && stm.getExpression() instanceof BinaryExpression ) {
            def binary = (BinaryExpression)stm.getExpression()
            if(!(binary.leftExpression instanceof VariableExpression))
                return
            if( binary.operation.type != Types.MINUS )
                return

            // -- transform the binary expression into a method call expression
            //    where the left expression represents the method name to invoke
            def methodName = ((VariableExpression)binary.leftExpression).name

            // -- wrap the value into a minus operator
            def value = (Expression)new UnaryMinusExpression( binary.rightExpression )
            def args = new ArgumentListExpression( [value] )

            // -- create the method call expression and replace it to the binary expression
            def call = new MethodCallExpression(new VariableExpression('this'), methodName, args)
            stm.setExpression(call)

        }
    }

    protected void fixStdinStdout( ExpressionStatement stm ) {

        // transform the following syntax:
        //      `stdin from x`  --> stdin() from (x)
        //      `stdout into x` --> `stdout() into (x)`
        if( stm.expression instanceof PropertyExpression ) {
            def expr = (PropertyExpression)stm.expression
            def obj = expr.objectExpression
            def prop = expr.property as ConstantExpression
            def target = new VariableExpression(prop.text)

            if( obj instanceof MethodCallExpression ) {
                def methodCall = obj as MethodCallExpression
                if( 'stdout' == methodCall.getMethodAsString() ) {
                    def stdout = new MethodCallExpression( new VariableExpression('this'), 'stdout', new ArgumentListExpression()  )
                    def into = new MethodCallExpression(stdout, 'into', new ArgumentListExpression(target))
                    // remove replace the old one with the new one
                    stm.setExpression( into )
                }
                else if( 'stdin' == methodCall.getMethodAsString() ) {
                    def stdin = new MethodCallExpression( new VariableExpression('this'), 'stdin', new ArgumentListExpression()  )
                    def from = new MethodCallExpression(stdin, 'from', new ArgumentListExpression(target))
                    // remove replace the old one with the new one
                    stm.setExpression( from )
                }
            }
        }
        // transform the following syntax:
        //      `stdout into (x,y,..)` --> `stdout() into (x,y,..)`
        else if( stm.expression instanceof MethodCallExpression ) {
            def methodCall = (MethodCallExpression)stm.expression
            if( 'stdout' == methodCall.getMethodAsString() ) {
                def args = methodCall.getArguments()
                if( args instanceof ArgumentListExpression && args.getExpressions() && args.getExpression(0) instanceof MethodCallExpression ) {
                    def methodCall2 = (MethodCallExpression)args.getExpression(0)
                    def args2 = methodCall2.getArguments()
                    if( args2 instanceof ArgumentListExpression && methodCall2.methodAsString == 'into') {
                        def vars = args2.getExpressions()
                        def stdout = new MethodCallExpression( new VariableExpression('this'), 'stdout', new ArgumentListExpression()  )
                        def into = new MethodCallExpression(stdout, 'into', new ArgumentListExpression(vars))
                        // remove replace the old one with the new one
                        stm.setExpression( into )
                    }
                }
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

            if( methodName in ['val','env','file','each','set','stdin'] ) {
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

    private boolean withinEachMethod

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
        withinEachMethod = name == '_in_each'

        try {
            if( isOutputValWithPropertyExpression(methodCall) )
                // transform an output value declaration such
                //   output: val( obj.foo )
                // to
                //   output: val({ obj.foo })
                wrapPropertyToClosure((ArgumentListExpression)methodCall.getArguments())
            else
                varToConst(methodCall.getArguments())

        } finally {
            withinSetMethod = false
            withinFileMethod = false
            withinEachMethod = false
        }
    }

    protected boolean isOutputValWithPropertyExpression(MethodCallExpression methodCall) {
        if( methodCall.methodAsString != '_out_val' )
            return false
        if( methodCall.getArguments() instanceof ArgumentListExpression ) {
            def args = (ArgumentListExpression)methodCall.getArguments()
            if( args.size() != 1 )
                return false

            return args.getExpression(0) instanceof PropertyExpression
        }

        return false
    }

    protected void wrapPropertyToClosure(ArgumentListExpression expr) {
        def args = expr as ArgumentListExpression
        def property = (PropertyExpression) args.getExpression(0)

        def closure = wrapPropertyToClosure(property)

        args.getExpressions().set(0, closure)
    }

    protected ClosureExpression wrapPropertyToClosure(PropertyExpression property)  {
        def block = new BlockStatement()
        block.addStatement( new ExpressionStatement(property) )

        def closure = new ClosureExpression( Parameter.EMPTY_ARRAY, block )
        closure.variableScope = new VariableScope(block.variableScope)

        return closure
    }


    protected Expression varToStr( Expression expr ) {
        if( expr instanceof VariableExpression ) {
            def name = ((VariableExpression) expr).getName()
            return newObj( TokenVar, new ConstantExpression(name) )
        }
        else if( expr instanceof PropertyExpression ) {
            // transform an output declaration such
            // output: set val( obj.foo )
            //  to
            // output: set val({ obj.foo })
            return wrapPropertyToClosure(expr)
        }

        if( expr instanceof TupleExpression )  {
            def i = 0
            def list = expr.getExpressions()
            for( Expression item : list ) {
                list[i++] = varToStr(item)
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

        if( expr instanceof MethodCallExpression ) {
            def methodCall = expr as MethodCallExpression

            /*
             * replace 'file' method call in the set definition, for example:
             *
             * input:
             *   set( file(fasta:'*.fa'), .. ) from q
             */
            if( methodCall.methodAsString == 'file' && (withinSetMethod || withinEachMethod) ) {
                def args = (TupleExpression) varToConst(methodCall.arguments)
                return newObj( TokenFileCall, args )
            }

            /*
             * input:
             *  set( env(VAR_NAME) ) from q
             */
            if( methodCall.methodAsString == 'env' && withinSetMethod ) {
                def args = (TupleExpression) varToStr(methodCall.arguments)
                return newObj( TokenEnvCall, args )
            }

            /*
             * input:
             *   set val(x), .. from q
             */
            if( methodCall.methodAsString == 'val' && withinSetMethod ) {
                def args = (TupleExpression) varToStr(methodCall.arguments)
                return newObj( TokenValCall, args )
            }

        }

        // -- TupleExpression or ArgumentListExpression
        if( expr instanceof TupleExpression )  {
            def i = 0
            def list = expr.getExpressions()
            for( Expression item : list )  {
                list[i++] = varToConst(item)
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
    protected boolean wrapExpressionWithClosure( BlockStatement block, Expression expr, int len, CharSequence source, SourceUnit unit ) {
        if( expr instanceof GStringExpression || expr instanceof ConstantExpression ) {
            // remove the last expression
            block.statements.remove(len-1)

            // and replace it by a wrapping closure
            def closureExp = new ClosureExpression( Parameter.EMPTY_ARRAY, new ExpressionStatement(expr) )
            closureExp.variableScope = new VariableScope(block.variableScope)

            // append to the list of statement
            //def wrap = newObj(TaskBody, closureExp, new ConstantExpression(source.toString()), ConstantExpression.TRUE)
            def wrap = makeScriptWrapper(closureExp, source, 'script', unit )
            block.statements.add( new ExpressionStatement(wrap) )

            return true
        }
        else if( expr instanceof ClosureExpression ) {
            // do not touch it
            return true
        }
        else {
            log.trace "Invalid process result expression: ${expr} -- Only constant or string expression can be used"
        }

        return false
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
        if( list.size() != 1 || !list[0].class.isAssignableFrom(MethodCallExpression) ) {
            log.debug "Missing name in process definition at line: ${methodCall.lineNumber}"
            unit.addError( new SyntaxException("Process definition syntax error -- You must provide a string identifier after the `process` keyword", methodCall.lineNumber, methodCall.columnNumber+7))
            return
        }

        def nested = list[0] as MethodCallExpression
        def name = nested.getMethodAsString()
        // check the process name is not defined yet
        if( !processNames.add(name) ) {
            log.warn "Process `$name` is defined two or more times"
        }

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
     * Fetch all the variable references in a closure expression.
     *
     * @param closure
     * @param unit
     * @return The set of variable names referenced in the script. NOTE: it includes properties in the form {@code object.propertyName}
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

        final Set<String> localDef = []

        final SourceUnit sourceUnit

        private boolean declaration

        VariableVisitor( SourceUnit unit ) {
            this.sourceUnit = unit
        }

        protected boolean isNormalized(PropertyExpression expr) {
            if( !(expr.getProperty() instanceof ConstantExpression) )
                return false

            def target = expr.getObjectExpression()
            while( target instanceof PropertyExpression) {
                target = (target as PropertyExpression).getObjectExpression()
            }

            return target instanceof VariableExpression
        }

        @Override
        void visitDeclarationExpression(DeclarationExpression expr) {
            declaration = true
            try {
                super.visitDeclarationExpression(expr)
            }
            finally {
                declaration = false
            }
        }

        @Override
        void visitPropertyExpression(PropertyExpression expr) {

            if( isNormalized(expr)) {
                final name = expr.text.replace('?','')
                final line = expr.lineNumber
                final coln = expr.columnNumber

                if( !name.startsWith('this.') && !fAllVariables.containsKey(name) ) {
                    fAllVariables[name] = new TokenValRef(name,line,coln)
                }
            }
            else
                super.visitPropertyExpression(expr)

        }

        @Override
        void visitVariableExpression(VariableExpression var) {
            final name = var.name
            final line = var.lineNumber
            final coln = var.columnNumber

            if( name == 'this' )
                return

            if( declaration ) {
                if( fAllVariables.containsKey(name) )
                    sourceUnit.addError( new SyntaxException("Variable `$name` already defined in the process scope", line, coln))
                else
                    localDef.add(name)
            }

            // Note: variable declared in the process scope are not added
            // to the set of referenced variables. Only global ones are tracked
            else if( !localDef.contains(name) && !fAllVariables.containsKey(name) ) {
                fAllVariables[name] = new TokenValRef(name,line,coln)
            }
        }

        @Override
        protected SourceUnit getSourceUnit() {
            return sourceUnit
        }

        /**
         * @return The set of all variables referenced in the script.
         * NOTE: it includes properties in the form {@code object.propertyName}
         */
        Set<TokenValRef> getAllVariables() {
            new HashSet<TokenValRef>(fAllVariables.values())
        }
    }

    /**
     * Transform any GString to a Lazy GString i.e.
     *
     * from
     *   "${foo} ${bar}"
     * to
     *   "${->foo} ${->bar}
     *
     */
    static class GStringToLazyVisitor extends ClassCodeVisitorSupport {

        final SourceUnit sourceUnit

        private withinClosure

        GStringToLazyVisitor(SourceUnit unit) {
            this.sourceUnit = unit
        }

        @Override
        public void visitClosureExpression(ClosureExpression expression) {
            withinClosure = true
            try {
                super.visitClosureExpression(expression)
            }
            finally {
                withinClosure = false
            }
        }

        @Override
        void visitGStringExpression(GStringExpression expression) {
            if( !withinClosure ) {
                xformToLazy(expression)
            }
        }

        protected void xformToLazy(GStringExpression str) {
            def values = str.getValues()
            def normalised = new Expression[values.size()]

            // wrap all non-closure to a ClosureExpression
            for( int i=0; i<values.size(); i++ ) {
                final item = values[i]
                if( item instanceof ClosureExpression  ) {
                    // when there is already a closure the conversion it is aborted
                    // because it supposed the gstring is already a lazy-string
                    return
                }
                normalised[i] = wrapWithClosure(item)
            }

            for( int i=0; i<values.size(); i++ ) {
                values[i] = normalised[i]
            }
        }

        protected ClosureExpression wrapWithClosure( Expression expr ) {

            // create an expression statement for the given `expression`
            def statement = new ExpressionStatement(expr)
            // add it to a new block
            def block = new BlockStatement()
            block.addStatement(statement)
            // create a closure over the given block
            // note: the closure parameter argument must be *null* to force the creation of a closure like {-> something}
            // otherwise it creates a closure with an implicit parameter that is managed in a different manner by the
            // GString -- see http://docs.groovy-lang.org/latest/html/documentation/#_special_case_of_interpolating_closure_expressions
            new ClosureExpression( null, block )

        }

        @Override
        protected SourceUnit getSourceUnit() {
            return sourceUnit
        }

    }

}
