/*
 * Copyright 2013-2023, Seqera Labs
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

import static nextflow.Const.*
import static nextflow.ast.ASTHelpers.*
import static org.codehaus.groovy.ast.tools.GeneralUtils.*

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.NF
import nextflow.script.BaseScript
import nextflow.script.BodyDef
import nextflow.script.IncludeDef
import nextflow.script.TaskClosure
import nextflow.script.TokenEnvCall
import nextflow.script.TokenFileCall
import nextflow.script.TokenPathCall
import nextflow.script.TokenStdinCall
import nextflow.script.TokenStdoutCall
import nextflow.script.TokenValCall
import nextflow.script.TokenValRef
import nextflow.script.TokenVar
import org.codehaus.groovy.ast.ASTNode
import org.codehaus.groovy.ast.ClassCodeVisitorSupport
import org.codehaus.groovy.ast.ClassNode
import org.codehaus.groovy.ast.ConstructorNode
import org.codehaus.groovy.ast.MethodNode
import org.codehaus.groovy.ast.Parameter
import org.codehaus.groovy.ast.VariableScope
import org.codehaus.groovy.ast.expr.ArgumentListExpression
import org.codehaus.groovy.ast.expr.BinaryExpression
import org.codehaus.groovy.ast.expr.CastExpression
import org.codehaus.groovy.ast.expr.ClosureExpression
import org.codehaus.groovy.ast.expr.ConstantExpression
import org.codehaus.groovy.ast.expr.Expression
import org.codehaus.groovy.ast.expr.GStringExpression
import org.codehaus.groovy.ast.expr.ListExpression
import org.codehaus.groovy.ast.expr.MapEntryExpression
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
import org.codehaus.groovy.syntax.Token
import org.codehaus.groovy.syntax.Types
import org.codehaus.groovy.transform.ASTTransformation
import org.codehaus.groovy.transform.GroovyASTTransformation

/**
 * Implement the AST transformations for the Nextflow DSL.
 *
 * Includes transforms for functions, imports, processes, and workflows.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
@GroovyASTTransformation(phase = CompilePhase.CONVERSION)
class NextflowDSLImpl implements ASTTransformation {

    final static private String WORKFLOW_TAKE = 'take'
    final static private String WORKFLOW_EMIT = 'emit'
    final static private String WORKFLOW_MAIN = 'main'
    final static private List<String> SCOPES = [WORKFLOW_TAKE, WORKFLOW_EMIT, WORKFLOW_MAIN]

    final static public String PROCESS_WHEN = 'when'
    final static public String PROCESS_STUB = 'stub'

    static public String OUT_PREFIX = '$out'

    static private Set<String> RESERVED_NAMES

    static {
        // method names implicitly defined by the groovy script SHELL
        RESERVED_NAMES = ['main','run','runScript'] as Set
        // existing method cannot be used for custom script definition
        for( def method : BaseScript.getMethods() ) {
            RESERVED_NAMES.add(method.name)
        }
    }

    @Override
    void visit(ASTNode[] astNodes, SourceUnit unit) {
        createVisitor(unit).visitClass((ClassNode)astNodes[1])
    }

    protected ClassCodeVisitorSupport createVisitor( SourceUnit unit ) {
        new DslCodeVisitor(unit)
    }

    @CompileStatic
    static class DslCodeVisitor extends ClassCodeVisitorSupport {

        final private SourceUnit unit

        private Set<String> processNames = []

        private Set<String> workflowNames = []

        private Set<String> functionNames = []

        private int anonymousWorkflow

        protected SourceUnit getSourceUnit() { unit }

        DslCodeVisitor(SourceUnit unit) {
            this.unit = unit
        }

        /**
         * Intercept method declarations to register them as function names.
         *
         * @param method
         */
        @Override
        void visitMethod(MethodNode method) {
            if( method.public && !method.static && !method.synthetic && !method.metaDataMap?.'org.codehaus.groovy.ast.MethodNode.isScriptBody' ) {
                if( !isIllegalName(method.name, method) )
                    functionNames.add(method.name)
            }
            super.visitMethod(method)
        }

        /**
         * Intercept top-level method calls to 'process' and 'workflow' and
         * transform them into process and workflow definitions.
         *
         * @param methodCall
         */
        @Override
        void visitMethodCallExpression(MethodCallExpression methodCall) {
            // determine whether the method call is in the top-level scope
            final isTopLevel = methodCall.objectExpression?.getText() == 'this'
            final methodName = methodCall.getMethodAsString()

            // transform 'process' method call into process definition
            if( methodName == 'process' && isTopLevel ) {
                convertProcessDef(methodCall,sourceUnit)
            }

            // transform 'workflow' method call into workflow definition
            else if( methodName == 'workflow' && isTopLevel ) {
                convertWorkflowDef(methodCall,sourceUnit)
            }

            super.visitMethodCallExpression(methodCall)
        }

        /**
         * Intercept top-level method calls to 'include' as module imports.
         *
         * For example, the following module import:
         *
         *   include { foo; foo as bar } from './some/module' addParams(foo: 'Ciao')
         *
         * is parsed by Groovy into:
         *
         *   this.include({ foo; foo as bar }).from('./some/module').addParams([foo: 'Ciao'])
         *
         * and then transformed into:
         *
         *   this.include( IncludeDef([ Module('foo'), Module('foo', 'bar') ]) )
         *       .from('./some/module')
         *       .addParams([foo: 'Ciao'])
         *       .load0(params)
         *
         * @param stm
         */
        @Override
        void visitExpressionStatement(ExpressionStatement stm) {
            if( stm.text.startsWith('this.include(') && stm.getExpression() instanceof MethodCallExpression )  {
                // transform method call into module import
                final methodCall = (MethodCallExpression)stm.getExpression()
                convertIncludeDef(methodCall)
                // invoke load0() (from IncludeDef) on method call expression
                final loadCall = new MethodCallExpression(methodCall, 'load0', new ArgumentListExpression(new VariableExpression('params')))
                stm.setExpression(loadCall)
            }
            super.visitExpressionStatement(stm)
        }

        /**
         * Transform an 'include' method call into a module import.
         *
         * @param call
         */
        protected void convertIncludeDef(MethodCallExpression call) {
            if( call.methodAsString=='include' && call.arguments instanceof ArgumentListExpression ) {
                // make sure there is a single closure argument
                final allArgs = (ArgumentListExpression)call.arguments
                if( allArgs.size() != 1 || allArgs[0] !instanceof ClosureExpression ) {
                    syntaxError(call, "Not a valid include definition -- it must specify a single closure argument")
                    return
                }

                // extract module arguments from closure
                final arg = (ClosureExpression)allArgs[0]
                final block = (BlockStatement)arg.getCode()
                final modulesList = new ListExpression()
                for( Statement stm : block.statements ) {
                    if( stm instanceof ExpressionStatement ) {
                        CastExpression castX
                        VariableExpression varX
                        Expression moduleX
                        // extract module name, e.g. `foo`
                        if( (varX=isVariableX(stm.expression)) ) {
                            def name = constX(varX.name)
                            moduleX = createX(IncludeDef.Module, name)
                        }
                        // extract module name with alias, e.g. `foo as bar`
                        else if( (castX=isCastX(stm.expression)) && (varX=isVariableX(castX.expression)) ) {
                            def name = constX(varX.name)
                            final alias = constX(castX.type.name)
                            moduleX = createX(IncludeDef.Module, name, alias)
                        }
                        // otherwise return an error
                        else {
                            syntaxError(call, "Not a valid include module name")
                            return
                        }
                        modulesList.addExpression(moduleX)
                    }
                    else {
                        syntaxError(call, "Not a valid include module name")
                        return
                    }
                }

                // replace include() argument with IncludeDef instance
                call.setArguments(new ArgumentListExpression( createX(IncludeDef, modulesList) ))
            }
            else if( call.objectExpression instanceof MethodCallExpression ) {
                convertIncludeDef((MethodCallExpression)call.objectExpression)
            }
        }

        /**
         * Transform a 'workflow' method call into a workflow definition.
         *
         * For example, the following Nextflow code:
         *
         *   workflow foo { code }
         *
         * is parsed by Groovy into:
         *
         *   this.workflow( foo({ code }) )
         *
         * and then transformed into:
         *
         *   this.workflow( 'foo', { new BodyDef({ code }, 'code') } )
         *
         * @param methodCall
         * @param unit
         */
        protected void convertWorkflowDef(MethodCallExpression methodCall, SourceUnit unit) {
            log.trace "Convert 'workflow' ${methodCall.arguments}"

            assert methodCall.arguments instanceof ArgumentListExpression
            def args = (ArgumentListExpression)methodCall.arguments

            // create anonymous workflow definition if there is a single closure argument
            if( args.size() == 1 && args[0] instanceof ClosureExpression ) {
                // make sure there is only one anonymous workflow
                if( anonymousWorkflow++ > 0 ) {
                    unit.addError( new SyntaxException("Duplicate anonymous workflow definition", methodCall.lineNumber, methodCall.columnNumber+8) )
                    return
                }

                // extract the workflow body
                def body = (ClosureExpression)args[0]
                methodCall.setArguments(
                    new ArgumentListExpression( makeWorkflowBodyDefWrapper(body,true) )
                )
                return 
            }

            // otherwise, make sure there is a single argument, which is a method call expression
            if( args.size() != 1 || !args[0].class.isAssignableFrom(MethodCallExpression) ) {
                log.debug "Missing name in workflow definition at line: ${methodCall.lineNumber}"
                unit.addError( new SyntaxException("Workflow definition syntax error -- A string identifier must be provided after the `workflow` keyword", methodCall.lineNumber, methodCall.columnNumber+8))
                return
            }

            // the nested method name is the workflow name
            final nested = args[0] as MethodCallExpression
            final name = nested.getMethodAsString()

            // make sure the name is not reserved or already defined
            if( isIllegalName(name, methodCall) ) {
                return
            }

            // register the workflow name
            workflowNames.add(name)

            // make sure there is a single nested argument, which is a closure
            args = (ArgumentListExpression)nested.getArguments()
            log.trace "Workflow name: $name with args: $args"

            if( args.size() != 1 || args[0] !instanceof ClosureExpression ) {
                syntaxError(methodCall, "Invalid workflow definition")
                return
            }

            // extract the workflow body
            final body = (ClosureExpression)args[0]
            methodCall.setArguments(
                new ArgumentListExpression( constX(name), makeWorkflowBodyDefWrapper(body,false) )
            )
        }

        /**
         * TODO: normalize workflow params (take and emit)
         *
         * @param stat
         * @param type
         * @param uniqueNames
         * @param body
         */
        protected Statement normWorkflowParam(ExpressionStatement stat, String type, Set<String> uniqueNames, List<Statement> body) {
            MethodCallExpression callx
            VariableExpression varx

            if( (callx=isMethodCallX(stat.expression)) && isThisX(callx.objectExpression) ) {
                final name = "_${type}_${callx.methodAsString}"
                return stmt( callThisX(name, callx.arguments) )
            }

            if( (varx=isVariableX(stat.expression)) ) {
                final name = "_${type}_${varx.name}"
                return stmt( callThisX(name) )
            }

            if( type == WORKFLOW_EMIT ) {
                return createAssignX(stat, body, type, uniqueNames)
            }

            syntaxError(stat, "Workflow malformed parameter definition")
            return stat
        }

        /**
         * TODO: create assignment expression
         *
         * @param stat
         * @param body
         * @param type
         * @param uniqueNames
         */
        protected Statement createAssignX(ExpressionStatement stat, List<Statement> body, String type, Set<String> uniqueNames) {
            BinaryExpression binx
            MethodCallExpression callx
            Expression args=null

            if( (binx=isAssignX(stat.expression)) ) {
                // keep the statement in body to allow it to be evaluated
                body.add(stat)
                // and create method call expr to capture the var name in the emission
                final left = (VariableExpression)binx.leftExpression
                final name = "_${type}_${left.name}"
                return stmt( callThisX(name) )
            }

            if( (callx=isMethodCallX(stat.expression)) && callx.objectExpression.text!='this' && hasTo(callx)) {
                // keep the args
                args = callx.arguments
                // replace the method call expression with a property
                stat.expression = new PropertyExpression(callx.objectExpression, callx.method)
                // then, fallback to default case
            }

            // wrap the expression into a assignment expression
            final var = getNextName(uniqueNames)
            final left = new VariableExpression(var)
            final right = stat.expression
            final token = new Token(Types.ASSIGN, '=', -1, -1)
            final assign = new BinaryExpression(left, token, right)
            body.add(stmt(assign))

            // the call method statement for the emit declaration
            final name="_${type}_${var}"
            callx =  args ? callThisX(name, args) : callThisX(name)
            return stmt(callx)
        }

        protected boolean hasTo(MethodCallExpression callX) {
            def tupleX = isTupleX(callX.arguments)
            if( !tupleX ) return false
            if( !tupleX.expressions ) return false
            def mapX = isMapX(tupleX.expressions[0])
            if( !mapX ) return false
            def entry = mapX.getMapEntryExpressions().find { isConstX(it.keyExpression).text=='to' }
            return entry != null
        }

        protected String getNextName(Set<String> allNames) {
            String result
            while( true ) {
                result = OUT_PREFIX + allNames.size()
                if( allNames.add(result) )
                    break
            }
            return result
        }

        /**
         * Create a BodyDef closure from a workflow body (closure expression).
         *
         * @param closure
         * @param anonymous
         */
        protected Expression makeWorkflowBodyDefWrapper( ClosureExpression closure, boolean anonymous ) {

            final codeBlock = (BlockStatement) closure.code
            final codeStms = codeBlock.statements
            final scope = codeBlock.variableScope

            final emitNames = new LinkedHashSet<String>(codeStms.size())
            final wrap = new ArrayList<Statement>(codeStms.size())
            final body = new ArrayList<Statement>(codeStms.size())

            String currLabel = null
            String prevLabel = null
            final visited = new HashMap<String,Boolean>(5)
            for( Statement stm : codeStms ) {
                // update the current and previous label
                prevLabel = currLabel
                currLabel = stm.statementLabel ?: currLabel

                // return an error if a label is repeated after a different label
                if( currLabel && currLabel != prevLabel ) {
                    if( visited[currLabel] && visited[prevLabel] ) {
                        syntaxError(stm, "Unexpected workflow label `${currLabel}:` here")
                        break
                    }
                }
                visited[currLabel] = true

                switch (currLabel) {
                    // add statements after 'take:' and 'emit:' to the workflow inputs and outputs
                    case WORKFLOW_TAKE:
                    case WORKFLOW_EMIT:
                        if( !(stm instanceof ExpressionStatement) ) {
                            syntaxError(stm, "Workflow malformed parameter definition")
                            break
                        }
                        wrap.add(normWorkflowParam(stm as ExpressionStatement, currLabel, emitNames, body))
                    break

                    // add statements after 'main:' to the workflow body
                    case WORKFLOW_MAIN:
                        body.add(stm)
                        break

                    // add unlabeled statements to the workflow body
                    default:
                        if( currLabel ) {
                            def opts = SCOPES.closest(currLabel)
                            def msg = "Unknown workflow label `${currLabel}:`"
                            if( opts ) msg += " -- Did you mean ${opts.collect{"'$it'"}.join(', ')}"
                            syntaxError(stm, msg)
                        }
                        body.add(stm)
                }
            }

            // read the closure source into a string
            final source = new StringBuilder()
            readSource(closure, source, unit, true)

            // return new closure expression with extracted body
            final bodyClosure = closureX(null, block(scope, body))
            final invokeBody = makeBodyDefWrapper(bodyClosure, source.toString(), 'workflow', unit)
            wrap.add( stmt(invokeBody) )

            closureX(null, block(scope, wrap))
        }

        /**
         * Append a syntax error at the given AST node.
         *
         * @param node
         * @param message
         */
        protected void syntaxError(ASTNode node, String message) {
            int line = node.lineNumber
            int coln = node.columnNumber
            unit.addError( new SyntaxException(message,line,coln) )
        }

        /**
         * Transform a process body into a {@code BodyDef} closure.
         *
         * @param body
         */
        protected void convertProcessBody( BlockStatement body ) {

            List<Statement> execStatements = []
            def execSource = new StringBuilder()

            List<Statement> whenStatements = []
            def whenSource = new StringBuilder()

            List<Statement> stubStatements = []
            def stubSource = new StringBuilder()

            def iterator = body.getStatements().iterator()
            String currLabel = null
            String bodyLabel = null

            // iterate through each statement in the body closure
            while( iterator.hasNext() ) {

                // get next statement
                Statement stm = iterator.next()

                // update the current label
                currLabel = stm.statementLabel ?: currLabel

                switch(currLabel) {
                    // transform statements after the 'input:' label into input parameters
                    case 'input':
                        if( stm instanceof ExpressionStatement ) {
                            fixLazyGString( stm )
                            fixStdinStdout( stm )
                            convertInputMethod( stm.getExpression() )
                        }
                        break

                    // transform statements after the 'output:' label into output parameters
                    case 'output':
                        if( stm instanceof ExpressionStatement ) {
                            fixLazyGString( stm )
                            fixStdinStdout( stm )
                            convertOutputMethod( stm.getExpression() )
                        }
                        break

                    // collect all statements after the 'exec:', 'script:', and 'shell:' labels
                    // and remove them from the process body
                    case 'exec':
                    case 'script':
                    case 'shell':
                        bodyLabel = currLabel
                        iterator.remove()
                        execStatements << stm
                        readSource(stm,execSource,unit)
                        break

                    // collect all statements after the 'stub:' label
                    // and remove them from the process body
                    case PROCESS_STUB:
                        iterator.remove()
                        stubStatements << stm
                        readSource(stm,stubSource,unit)
                        break

                    // collect all statements after a 'when:' label
                    // and remove them from the process body
                    case PROCESS_WHEN:
                        if( iterator.hasNext() ) {
                            iterator.remove()
                            whenStatements << stm
                            readSource(stm,whenSource,unit)
                            break
                        }

                        // return a syntax error if a when statement is the last statement
                        // in the process body (should be the task command)
                        else if( !whenStatements ) {
                            int line = body.lineNumber
                            int coln = body.columnNumber
                            unit.addError(new SyntaxException("Invalid process definition -- Empty `when` or missing `script` statement", line, coln))
                            return
                        }

                        // otherwise, handle the next statement as the task command
                        else
                            break

                    // leave all other statements in the process body
                    default:
                        if( currLabel ) {
                            def line = stm.getLineNumber()
                            def coln = stm.getColumnNumber()
                            unit.addError(new SyntaxException("Invalid process definition -- Unknown label `$currLabel`",line,coln))
                            return
                        }
                        fixLazyGString(stm)
                        fixDirectiveWithNegativeValue(stm)
                }
            }

            // add the `when` block if found
            if( whenStatements ) {
                addNamedBlock(PROCESS_WHEN, whenStatements, whenSource, body)
            }

            // add the `stub` block if found
            if( stubStatements ) {
                final stubBlock = addNamedBlock(PROCESS_STUB, stubStatements, stubSource, body)
                stubBlock.visit(new TaskCmdXformVisitor(unit))
            }

            // attempt to create a BodyDef wrapper with the task command
            boolean done = false
            final len = body.statements.size()

            // use the `exec` block if found
            if( execStatements ) {
                // wrap the exec statements in a BodyDef closure
                def execClosure = closureX(
                    null,
                    block(body.variableScope, execStatements)
                )
                final wrap = makeBodyDefWrapper(execClosure, execSource, bodyLabel, unit)

                // append the BodyDef closure to the process body
                body.addStatement( new ExpressionStatement(wrap) )
                // set the 'script' flag parameter
                if( bodyLabel == 'script' )
                    body.visit(new TaskCmdXformVisitor(unit))
                done = true
            }

            // otherwise, add an empty command when only the `stub` block is defined
            else if ( !bodyLabel && stubStatements ) {
                // wrap an empty statement in a BodyDef closure
                final cmd = 'true'
                final dummyClosure = closureX(
                    null,
                    block(body.variableScope, [ new ExpressionStatement(constX(cmd)) ] as List<Statement>)
                )
                final wrap = makeBodyDefWrapper(dummyClosure, cmd, 'script', unit)

                // append the new block to the process body
                body.addStatement( new ExpressionStatement(wrap) )
                done = true
            }

            // otherwise, add the last statement as the task command if it exists
            // (it can be specified without the `script` label)
            else if( len ) {
                // get the last statement
                def stm = body.getStatements().get(len-1)
                readSource(stm,execSource,unit)

                // attempt to wrap the last statement in a BodyDef closure
                if ( stm instanceof ReturnStatement ) {
                    done = wrapLastStatementWithBodyDef(body, stm.getExpression(), len, execSource, unit)
                }

                else if ( stm instanceof ExpressionStatement ) {
                    done = wrapLastStatementWithBodyDef(body, stm.getExpression(), len, execSource, unit)
                }

                // apply command variables escape
                stm.visit(new TaskCmdXformVisitor(unit))
            }

            // return a syntax error if the task closure could not be created
            if ( !done ) {
                log.trace "Invalid 'process' definition -- Process must terminate with string expression"
                int line = body.lineNumber
                int coln = body.columnNumber
                unit.addError( new SyntaxException("Invalid process definition -- Make sure the process ends with a script wrapped by quote characters",line,coln) )
            }
        }

        /**
         * Transform a named block (e.g. `when`, `stub`) into a method call on {@code ProcessConfig}
         * and append it to the parent block.
         *
         * For example, the following code:
         *
         *   when: params.foo
         *   stub: 'true'
         *
         * is transformed into:
         *
         *   this.when( new TaskClosure({ params.foo }, 'params.foo') )
         *   this.stub( new TaskClosure({ 'true' }, 'true') )
         *
         * See {@link nextflow.script.ProcessConfig#configProperties}
         * See {@link nextflow.processor.TaskConfig#getGuard(java.lang.String)}
         *
         * @param blockName
         * @param statements
         * @param source
         * @param parent
         */
        protected BlockStatement addNamedBlock( String blockName, List<Statement> statements, StringBuilder source, BlockStatement parent ) {
            // wrap the statements in a closure expression
            def block = new BlockStatement(statements, new VariableScope(parent.variableScope))
            def closure = closureX( null, block )

            // wrap the closure in a {@code TaskClosure} in order to capture the closure source
            def taskClosure = createX( TaskClosure, [
                closure,
                new ConstantExpression(source.toString())
            ] )

            // creates a method call expression for the method `when`
            def methodCall = new MethodCallExpression(VariableExpression.THIS_EXPRESSION, blockName, taskClosure)

            // append the method call to the parent block
            parent.getStatements().add(0, new ExpressionStatement(methodCall))

            return block
        }

        /**
         * Create a {@code BodyDef} expression from a closure expression and source string.
         *
         * @param closure
         * @param source
         * @param section
         * @param unit
         */
        private Expression makeBodyDefWrapper( ClosureExpression closure, CharSequence source, String section, SourceUnit unit ) {

            // collect all variable tokens into a single list argument
            final variables = fetchVariables(closure,unit)
            final listArg = new ArrayList(variables.size())
            for( TokenValRef var : variables ) {
                listArg << createX( TokenValRef, [
                    new ConstantExpression(var.name),
                    new ConstantExpression(var.lineNum),
                    new ConstantExpression(var.colNum)
                ] as List<Expression> )
            }

            // create the BodyDef constructor call
            createX( BodyDef, [
                closure,
                new ConstantExpression(source.toString()),
                new ConstantExpression(section),
                new ListExpression(listArg)
            ] )
        }

        /**
         * Read the source code of an AST node into a string buffer.
         *
         * @param node
         * @param buffer
         * @param unit
         * @param stripBrackets
         */
        private void readSource( ASTNode node, StringBuilder buffer, SourceUnit unit, boolean stripBrackets=false ) {
            final colx = node.getColumnNumber()
            final colz = node.getLastColumnNumber()
            final first = node.getLineNumber()
            final last = node.getLastLineNumber()
            for( int i=first; i<=last; i++ ) {
                def line = unit.source.getLine(i, null)
                if( i==last ) {
                    line = line.substring(0,colz-1)
                    if( stripBrackets ) {
                        line = line.replaceFirst(/}.*$/,'')
                        if( !line.trim() ) continue
                    }
                }
                if( i==first ) {
                    line = line.substring(colx-1)
                    if( stripBrackets ) {
                        line = line.replaceFirst(/^.*\{/,'').trim()
                        if( !line.trim() ) continue
                    }
                }
                buffer.append(line) .append('\n')
            }
        }

        /**
         * Transform any GStrings in a method call into lazy GStrings.
         *
         * @param stm
         */
        protected void fixLazyGString( Statement stm ) {
            if( stm instanceof ExpressionStatement && stm.getExpression() instanceof MethodCallExpression ) {
                new GStringToLazyVisitor(unit).visitExpressionStatement(stm)
            }
        }

        /**
         * Correct process directives with negative values by converting them
         * from binary expressions (as they are so parsed by Groovy) to method calls.
         *
         * For example, the following Nextflow code:
         *
         *   maxErrors -1
         *
         * is parsed by Groovy as a binary expression:
         *
         *   maxErrors - 1
         *
         * and then corrected to:
         *
         *   maxErrors(-1)
         *
         * Fixes #180
         *
         * @param stm
         */
        protected void fixDirectiveWithNegativeValue( Statement stm ) {
            if( stm instanceof ExpressionStatement && stm.getExpression() instanceof BinaryExpression ) {
                def binary = (BinaryExpression)stm.getExpression()
                if( binary.leftExpression !instanceof VariableExpression )
                    return
                if( binary.operation.type != Types.MINUS )
                    return

                // -- extract the method name from the left-hand value
                def methodName = ((VariableExpression)binary.leftExpression).name

                // -- wrap the right-hand value in a minus expression
                def value = (Expression)new UnaryMinusExpression( binary.rightExpression )

                // -- replace the binary expression with a method call expression
                def call = new MethodCallExpression(
                    new VariableExpression('this'),
                    methodName,
                    new ArgumentListExpression( value )
                )

                stm.setExpression(call)
            }
        }

        /**
         * Transform process `stdin` and `stdout` declarations into method calls
         * on {@code ProcessConfig}.
         *
         * For example, the following Nextflow code:
         *
         *   input: stdin
         *   output: stdout
         *
         * is transformed into:
         *
         *   this._in_stdin()
         *   this._out_stdout()
         *
         * @param stm
         */
        protected void fixStdinStdout( ExpressionStatement stm ) {
            VariableExpression varX = isVariableX(stm.expression)
            if( varX && (varX.name == 'stdin' || varX.name == 'stdout') ) {
                final name = varX.name == 'stdin' ? '_in_stdin' : '_out_stdout'
                final call = new MethodCallExpression(
                    new VariableExpression('this'),
                    name,
                    new ArgumentListExpression()
                )
                stm.setExpression(call)
            }
        }

        /**
         * Transform process input declarations into method calls on {@code ProcessConfig}.
         *
         * For example, the following Nextflow code:
         *
         *   val v
         *   path 'p'
         *   tuple val(tv), path(tp)
         *   stdin
         *
         * is parsed by Groovy as:
         *
         *   val(v)
         *   path('p')
         *   tuple(val(tv), path(tp))
         *   stdin
         *
         * and then transformed into:
         *
         *   this._in_val( new TokenVar('v') )
         *   this._in_path( new TokenVar('p') )
         *   this._in_tuple( new TokenValCall(new TokenVar('tv')), new TokenPathCall(new TokenVar('tp')) )
         *   this._in_stdin()
         *
         * @param expression
         */
        protected void convertInputMethod( Expression expression ) {
            log.trace "convert > input expression: $expression"

            if( expression instanceof MethodCallExpression ) {
                def methodCall = (MethodCallExpression) expression
                def methodName = methodCall.getMethodAsString()
                def nested = methodCall.objectExpression instanceof MethodCallExpression

                log.trace "convert > input method: $methodName"

                if( methodName in ['val','env','file','each','set','stdin','path','tuple'] ) {
                    // convert type qualifier to the corresponding {@code ProcessConfig} method name
                    if( !nested )
                        methodCall.setMethod( new ConstantExpression('_in_' + methodName) )

                    // wrap variable declarations in corresponding constructor calls
                    fixMethodCall(methodCall)
                }

                // traverse nested method calls
                if( expression.objectExpression instanceof MethodCallExpression ) {
                    convertInputMethod(methodCall.objectExpression)
                }
            }

            else if( expression instanceof PropertyExpression ) {
                // traverse nested method calls
                if( expression.objectExpression instanceof MethodCallExpression ) {
                    convertInputMethod(expression.objectExpression)
                }
            }
        }

        /**
         * Wrap the `emit` option of an output declaration in a string literal.
         *
         * For example, the following code:
         *
         *   output: path 'foo', emit: bar
         *
         * is transformed into:
         *
         *   output: path 'foo', emit: 'bar'
         *
         * @param call
         */
        protected void fixOutEmitOption(MethodCallExpression call) {
            List<Expression> args = isTupleX(call.arguments)?.expressions
            if( !args ) return
            if( args.size()<2 && (args.size()!=1 || call.methodAsString!='_out_stdout')) return
            MapExpression map = isMapX(args[0])
            if( !map ) return
            for( int i=0; i<map.mapEntryExpressions.size(); i++ ) {
                final entry = map.mapEntryExpressions[i]
                final key = isConstX(entry.keyExpression)
                final val = isVariableX(entry.valueExpression)
                if( key?.text == 'emit' && val ) {
                    map.mapEntryExpressions[i] = new MapEntryExpression(key, constX(val.text))
                }
            }
        }

        /**
         * Transform process output declarations into method calls on {@code ProcessConfig}.
         *
         * For example, the following Nextflow code:
         *
         *   val v
         *   path 'p'
         *   tuple val(tv), path(tp)
         *   stdout
         *
         * is parsed by Groovy as:
         *
         *   val(v)
         *   path('p')
         *   tuple(val(tv), path(tp))
         *   stdout
         *
         * and then transformed into:
         *
         *   this._out_val( new TokenVar('v') )
         *   this._out_path( new TokenVar('p') )
         *   this._out_tuple( new TokenValCall(new TokenVar('tv')), new TokenPathCall(new TokenVar('tp')) )
         *   this._out_stdout()
         *
         * @param expression
         */
        protected void convertOutputMethod( Expression expression ) {
            log.trace "convert > output expression: $expression"

            if( expression !instanceof MethodCallExpression ) {
                return
            }

            def methodCall = (MethodCallExpression) expression
            def methodName = methodCall.getMethodAsString()
            def nested = methodCall.objectExpression instanceof MethodCallExpression

            log.trace "convert > output method: $methodName"

            if( methodName in ['val','env','file','set','stdout','path','tuple'] && !nested ) {
                // convert type qualifier to the corresponding {@code ProcessConfig} method name
                methodCall.setMethod( new ConstantExpression('_out_' + methodName) )

                // wrap variable declarations in corresponding constructor calls
                fixMethodCall(methodCall)

                // wrap emit value in string literal
                fixOutEmitOption(methodCall)
            }

            // traverse nested method calls
            if( methodCall.objectExpression instanceof MethodCallExpression ) {
                convertOutputMethod(methodCall.objectExpression)
            }
        }

        private boolean withinTupleMethod

        private boolean withinEachMethod

        /**
         * Transform a process input/output declaration by wrapping the name
         * in a {@code TokenVar} or a closure, in order to make it possible to
         * reference variables before they exist.
         *
         * For example, the following code:
         *
         *   this._in_val( foo )
         *   this._out_val( bar )
         *   this._out_val( obj.prop )
         *
         * is transformed into:
         *
         *   this._in_val( new TokenVar('foo') )
         *   this._out_val( new TokenVar('bar') )
         *   this._out_val({ obj.prop })
         *
         * @param methodCall
         */
        protected void fixMethodCall( MethodCallExpression methodCall ) {
            final name = methodCall.methodAsString

            withinTupleMethod = name == '_in_set' || name == '_out_set' || name == '_in_tuple' || name == '_out_tuple'
            withinEachMethod = name == '_in_each'

            try {
                // wrap the name in a closure (if it is an output with a property expression)
                if( isOutputWithPropertyExpression(methodCall) )
                    wrapPropertyToClosure((ArgumentListExpression)methodCall.getArguments())

                // otherwise wrap the name in a TokenVar
                else
                    varToConstX(methodCall.getArguments())
            }
            finally {
                withinTupleMethod = false
                withinEachMethod = false
            }
        }

        static final private List<String> OUT_PROPERTY_VALID_TYPES = ['_out_val', '_out_env', '_out_file', '_out_path']

        /**
         * Determine whether a method call is a process output declaration
         * with a property expression.
         *
         * @param methodCall
         */
        protected boolean isOutputWithPropertyExpression(MethodCallExpression methodCall) {
            if( methodCall.methodAsString !in OUT_PROPERTY_VALID_TYPES )
                return false

            if( methodCall.getArguments() !instanceof ArgumentListExpression )
                return false

            def args = (ArgumentListExpression)methodCall.getArguments()

            return args.size() > 0
                && args.size() <= 2
                && args.last() instanceof PropertyExpression
        }

        /**
         * Wrap the last argument in an argument list in a closure.
         *
         * @param args
         */
        protected void wrapPropertyToClosure(ArgumentListExpression args) {
            final property = (PropertyExpression) args.last()
            final closure = wrapPropertyToClosure(property)
            args.getExpressions().set(args.size()-1, closure)
        }

        /**
         * Wrap a property expression in a closure.
         *
         * @param property
         */
        protected ClosureExpression wrapPropertyToClosure(PropertyExpression property)  {
            def block = new BlockStatement()
            block.addStatement( new ExpressionStatement(property) )

            def closure = closureX( null, block )
            closure.variableScope = new VariableScope(block.variableScope)

            return closure
        }

        /**
         * Wrap an expression in an appropriate constructor call.
         *
         * For example, the following code:
         *
         *   val v
         *   tuple val( obj.foo )
         *
         * is transformed into:
         *
         *   val( new TokenVar('v') )
         *   tuple val({ obj.foo })
         *
         * @param expr
         */
        protected Expression varToStrX( Expression expr ) {
            // wrap variable expressions with {@code TokenVar}
            if( expr instanceof VariableExpression ) {
                def name = ((VariableExpression) expr).getName()
                return createX( TokenVar, new ConstantExpression(name) )
            }

            // wrap property expressions in a closure
            if( expr instanceof PropertyExpression ) {
                return wrapPropertyToClosure(expr)
            }

            // recursively traverse tuple expressions
            if( expr instanceof TupleExpression ) {
                def i = 0
                def args = expr.getExpressions()
                for( Expression item : args ) {
                    args[i++] = varToStrX(item)
                }

                return expr
            }

            // otherwise, return the input
            return expr
        }

        /**
         * Wrap a variable declaration in the appopriate constructor call.
         *
         * For example, the following code:
         *
         *   tuple( path(p), stdin )
         *
         * is transformed into:
         *
         *   tuple( new TokenPathCall(new TokenVar('p')), new TokenStdinCall() )
         *
         * @param expr
         */
        protected Expression varToConstX( Expression expr ) {

            if( expr instanceof VariableExpression ) {
                def name = ((VariableExpression) expr).getName()

                // for `stdin` within a `tuple`, wrap with {@code TokenStdinCall}
                if( name == 'stdin' && withinTupleMethod )
                    return createX( TokenStdinCall )

                // for `stdout` within a `tuple`, wrap with {@code TokenStdoutCall}
                else if ( name == 'stdout' && withinTupleMethod )
                    return createX( TokenStdoutCall )

                // otherwise, wrap with {@code TokenVar}
                else
                    return createX( TokenVar, new ConstantExpression(name) )
            }

            if( expr instanceof MethodCallExpression ) {
                def methodCall = expr as MethodCallExpression

                // for `file` within a `tuple` or `each`, wrap with {@code TokenFileCall}
                if( methodCall.methodAsString == 'file' && (withinTupleMethod || withinEachMethod) ) {
                    def args = (TupleExpression) varToConstX(methodCall.arguments)
                    return createX( TokenFileCall, args )
                }

                // for `path` within a `tuple` or `each`, wrap with {@code TokenPathCall}
                if( methodCall.methodAsString == 'path' && (withinTupleMethod || withinEachMethod) ) {
                    def args = (TupleExpression) varToConstX(methodCall.arguments)
                    return createX( TokenPathCall, args )
                }

                // for `env` within a `tuple`, wrap with {@code TokenEnvCall}
                if( methodCall.methodAsString == 'env' && withinTupleMethod ) {
                    def args = (TupleExpression) varToStrX(methodCall.arguments)
                    return createX( TokenEnvCall, args )
                }

                // for `val` within a `tuple`, wrap with {@code TokenValCall}
                if( methodCall.methodAsString == 'val' && withinTupleMethod ) {
                    def args = (TupleExpression) varToStrX(methodCall.arguments)
                    return createX( TokenValCall, args )
                }
            }

            // recursively traverse tuple or argument list expressions
            if( expr instanceof TupleExpression )  {
                def i = 0
                def args = expr.getExpressions()
                for( Expression item : args )  {
                    args[i++] = varToConstX(item)
                }
                return expr
            }

            return expr
        }

        /**
         * Wrap a statment with a BodyDef and append it to a block. The statement must be
         * a GString or constant expression, or it must already be a closure. Returns true
         * if the statement is valid.
         *
         * @param block The block to which the resulting closure has to be appended
         * @param expr The expression to the wrapped in a closure
         * @param len
         * @param source
         * @param unit
         */
        protected boolean wrapLastStatementWithBodyDef( BlockStatement block, Expression expr, int len, CharSequence source, SourceUnit unit ) {
            // wrap statement in BodyDef if it is a GString or constant
            if( expr instanceof GStringExpression || expr instanceof ConstantExpression ) {
                // remove the last statement
                block.statements.remove(len-1)

                // and replace it by a wrapping closure
                def closureExp = new ClosureExpression( Parameter.EMPTY_ARRAY, new ExpressionStatement(expr) )
                closureExp.variableScope = new VariableScope(block.variableScope)

                // append to the block statement
                def wrap = makeBodyDefWrapper(closureExp, source, 'script', unit)
                block.statements.add( new ExpressionStatement(wrap) )

                return true
            }

            // do nothing if the statement is already a closure
            else if( expr instanceof ClosureExpression ) {
                return true
            }

            log.trace "Invalid process result expression: ${expr} -- Only constant or string expression can be used"
            return false
        }

        /**
         * Determine whether a name is reserved, or already defined as a
         * workflow or process, or contains a colon ':' character.
         *
         * @param name
         * @param node
         */
        protected boolean isIllegalName(String name, ASTNode node) {
            if( name in RESERVED_NAMES ) {
                unit.addError( new SyntaxException("Identifier `$name` is reserved for internal use", node.lineNumber, node.columnNumber + 8) )
                return true
            }
            if( name in workflowNames || name in processNames ) {
                unit.addError( new SyntaxException("Identifier `$name` is already used by another definition", node.lineNumber, node.columnNumber + 8) )
                return true
            }
            if( name.contains(SCOPE_SEP) ) {
                def offset = 8 + 2 + name.indexOf(SCOPE_SEP)
                unit.addError( new SyntaxException("Process and workflow names cannot contain colon character", node.lineNumber, node.columnNumber + offset) )
                return true
            }
            return false
        }

        /**
         * Transform a 'process' method call into a process definition.
         *
         * For example, the following Nextflow code:
         *
         *   process foo { code }
         *
         * is parsed by Groovy into:
         *
         *   this.process( foo({ code }) )
         *
         * and then transformed into:
         *
         *   this.process( 'foo', { new BodyDef({ code }, 'code') } )
         *
         * @param methodCall
         * @param unit
         */
        protected void convertProcessDef( MethodCallExpression methodCall, SourceUnit unit ) {
            log.trace "Convert 'process' ${methodCall.arguments}"

            assert methodCall.arguments instanceof ArgumentListExpression
            def args = (methodCall.arguments as ArgumentListExpression).getExpressions()

            // make sure there is a single argument, which is a method call expression
            if( args.size() != 1 || !args[0].class.isAssignableFrom(MethodCallExpression) ) {
                log.debug "Missing name in process definition at line: ${methodCall.lineNumber}"
                unit.addError( new SyntaxException("Process definition syntax error -- A string identifier must be provided after the `process` keyword", methodCall.lineNumber, methodCall.columnNumber+7) )
                return
            }

            // the nested method name is the process name
            def nested = args[0] as MethodCallExpression
            def name = nested.getMethodAsString()

            // make sure the name is not reserved or already defined
            if( isIllegalName(name, methodCall) ) {
                return
            }

            // register the process name
            processNames.add(name)

            // make sure there is a single nested argument, which is a closure
            def nestedArgs = (ArgumentListExpression)nested.getArguments()
            log.trace "Process name: $name with args: $nestedArgs"

            if( nestedArgs.size() != 1 || nestedArgs[0] !instanceof ClosureExpression ) {
                syntaxError(methodCall, "Invalid process definition")
                return
            }

            // transform the process body
            final body = (BlockStatement) ((ClosureExpression) nestedArgs[0]).code
            convertProcessBody(body)
        }

        /**
         * Fetch all variable references in a closure expression.
         *
         * NOTE: includes property expressions, e.g. `object.propertyName`
         *
         * @param closure
         * @param unit
         */
        protected Set<TokenValRef> fetchVariables( ClosureExpression closure, SourceUnit unit ) {
            def visitor = new VariableVisitor(unit)
            visitor.visitClosureExpression(closure)
            return visitor.allVariables
        }

    }

}
