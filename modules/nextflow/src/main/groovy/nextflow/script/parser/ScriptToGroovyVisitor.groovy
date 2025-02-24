/*
 * Copyright 2024, Seqera Labs
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
package nextflow.script.parser

import java.util.stream.Collectors

import groovy.transform.CompileStatic
import nextflow.ast.GStringToLazyVisitor
import nextflow.ast.TaskCmdXformVisitor
import nextflow.script.BaseScript
import nextflow.script.BodyDef
import nextflow.script.IncludeDef
import nextflow.script.TaskClosure
import nextflow.script.TokenEnvCall
import nextflow.script.TokenEvalCall
import nextflow.script.TokenFileCall
import nextflow.script.TokenPathCall
import nextflow.script.TokenStdinCall
import nextflow.script.TokenStdoutCall
import nextflow.script.TokenValCall
import nextflow.script.TokenVar
import nextflow.script.ast.AssignmentExpression
import nextflow.script.ast.FeatureFlagNode
import nextflow.script.ast.FunctionNode
import nextflow.script.ast.IncludeNode
import nextflow.script.ast.OutputNode
import nextflow.script.ast.ParamNode
import nextflow.script.ast.ProcessNode
import nextflow.script.ast.ScriptNode
import nextflow.script.ast.ScriptVisitorSupport
import nextflow.script.ast.WorkflowNode
import org.codehaus.groovy.ast.ASTNode
import org.codehaus.groovy.ast.VariableScope
import org.codehaus.groovy.ast.expr.ArgumentListExpression
import org.codehaus.groovy.ast.expr.BinaryExpression
import org.codehaus.groovy.ast.expr.ClosureExpression
import org.codehaus.groovy.ast.expr.DeclarationExpression
import org.codehaus.groovy.ast.expr.EmptyExpression
import org.codehaus.groovy.ast.expr.Expression
import org.codehaus.groovy.ast.expr.MethodCallExpression
import org.codehaus.groovy.ast.expr.PropertyExpression
import org.codehaus.groovy.ast.expr.TupleExpression
import org.codehaus.groovy.ast.expr.VariableExpression
import org.codehaus.groovy.ast.stmt.BlockStatement
import org.codehaus.groovy.ast.stmt.EmptyStatement
import org.codehaus.groovy.ast.stmt.ExpressionStatement
import org.codehaus.groovy.ast.stmt.ReturnStatement
import org.codehaus.groovy.ast.stmt.Statement
import org.codehaus.groovy.control.SourceUnit
import org.codehaus.groovy.syntax.SyntaxException

import static nextflow.script.ast.ASTHelpers.*
import static org.codehaus.groovy.ast.tools.GeneralUtils.*

/**
 * Visitor to convert a Nextflow script AST into a
 * Groovy AST which is executed against {@link BaseScript}.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
public class ScriptToGroovyVisitor extends ScriptVisitorSupport {

    private static Set<String> RESERVED_NAMES

    static {
        // internal groovy shell methods
        RESERVED_NAMES = ["main", "run", "runScript"] as Set
        // internal script dsl methods
        for( final method : BaseScript.class.getMethods() ) {
            RESERVED_NAMES.add(method.getName())
        }
    }

    private SourceUnit sourceUnit

    private ScriptNode moduleNode

    public ScriptToGroovyVisitor(SourceUnit sourceUnit) {
        this.sourceUnit = sourceUnit
        this.moduleNode = (ScriptNode) sourceUnit.getAST()
    }

    @Override
    protected SourceUnit getSourceUnit() {
        return sourceUnit
    }

    public void visit() {
        if( moduleNode == null )
            return
        super.visit(moduleNode)
        if( moduleNode.isEmpty() )
            moduleNode.addStatement(ReturnStatement.RETURN_NULL_OR_VOID)
    }

    @Override
    public void visitFeatureFlag(FeatureFlagNode node) {
        final names = node.name.tokenize('.')
        Expression target = varX(names.head())
        for( final name : names.tail() )
            target = propX(target, name)

        final result = stmt(assignX(target, node.value))
        moduleNode.addStatement(result)
    }

    @Override
    public void visitInclude(IncludeNode node) {
        final moduleArgs = (List<Expression>) node.modules.stream()
            .map((module) -> {
                final name = constX(module.name)
                return module.alias != null
                    ? createX(IncludeDef.Module.class, name, constX(module.alias))
                    : createX(IncludeDef.Module.class, name)
            })
            .collect(Collectors.toList())

        final include = callThisX("include", args(createX(IncludeDef.class, (Expression) listX(moduleArgs))))
        final from = callX(include, "from", args(node.source))
        final result = stmt(callX(from, "load0", args(varX("params"))))
        moduleNode.addStatement(result)
    }

    @Override
    public void visitParam(ParamNode node) {
        final result = stmt(assignX(node.target, node.value))
        moduleNode.addStatement(result)
    }

    @Override
    public void visitWorkflow(WorkflowNode node) {
        visitWorkflowTakes(node.takes)
        visit(node.main)
        visitWorkflowEmits(node.emits, node.main)
        visitWorkflowPublishers(node.publishers, node.main)

        final bodyDef = stmt(createX(
            BodyDef.class,
            args(
                closureX(node.main),
                constX(getSourceText(node)),
                constX("workflow")
            )
        ))
        final closure = closureX(block(new VariableScope(), List.of(
            node.takes,
            node.emits,
            bodyDef
        )))
        final arguments = node.isEntry()
            ? args(closure)
            : args(constX(node.getName()), closure)
        final result = stmt(callThisX("workflow", arguments))
        moduleNode.addStatement(result)
    }

    private void visitWorkflowTakes(Statement takes) {
        for( final stmt : asBlockStatements(takes) ) {
            final stmtX = (ExpressionStatement)stmt
            final take = (VariableExpression)stmtX.getExpression()
            stmtX.setExpression(callThisX("_take_", args(constX(take.getName()))))
        }
    }

    private void visitWorkflowEmits(Statement emits, Statement main) {
        final code = (BlockStatement)main
        for( final stmt : asBlockStatements(emits) ) {
            final stmtX = (ExpressionStatement)stmt
            final emit = stmtX.getExpression()
            if( emit instanceof VariableExpression ) {
                stmtX.setExpression(callThisX("_emit_", args(constX(emit.getName()))))
            }
            else if( emit instanceof AssignmentExpression ) {
                final target = (VariableExpression)emit.getLeftExpression()
                stmtX.setExpression(callThisX("_emit_", args(constX(target.getName()))))
                code.addStatement(stmtX)
            }
            else {
                final target = varX('$out')
                code.addStatement(assignS(target, emit))
                stmtX.setExpression(callThisX("_emit_", args(constX(target.getName()))))
                code.addStatement(stmtX)
            }
        }
    }

    private void visitWorkflowPublishers(Statement publishers, Statement main) {
        final code = (BlockStatement)main
        for( final stmt : asBlockStatements(publishers) ) {
            final stmtX = (ExpressionStatement)stmt
            final publish = (BinaryExpression)stmtX.getExpression()
            stmtX.setExpression(callThisX("_publish_target", args(publish.getLeftExpression(), publish.getRightExpression())))
            code.addStatement(stmtX)
        }
    }

    @Override
    public void visitProcess(ProcessNode node) {
        visitProcessDirectives(node.directives)
        visitProcessInputs(node.inputs)
        visitProcessOutputs(node.outputs)
        visit(node.exec)
        visit(node.stub)

        if( "script".equals(node.type) )
            node.exec.visit(new TaskCmdXformVisitor(sourceUnit))
        node.stub.visit(new TaskCmdXformVisitor(sourceUnit))

        final when = processWhen(node.when)
        final bodyDef = stmt(createX(
            BodyDef.class,
            args(
                closureX(node.exec),
                constX(getSourceText(node.exec)),
                constX(node.type)
            )
        ))
        final stub = processStub(node.stub)
        final closure = closureX(block(new VariableScope(), List.of(
            node.directives,
            node.inputs,
            node.outputs,
            when,
            stub,
            bodyDef
        )))
        final result = stmt(callThisX("process", args(constX(node.getName()), closure)))
        moduleNode.addStatement(result)
    }

    private void visitProcessDirectives(Statement directives) {
        asDirectives(directives).forEach((call) -> {
            fixLazyGString(call)
        })
    }

    private void visitProcessInputs(Statement inputs) {
        asDirectives(inputs).forEach((call) -> {
            fixLazyGString(call)

            final name = call.getMethodAsString()
            varToConstX(call.getArguments(), "tuple".equals(name), "each".equals(name))
            call.setMethod( constX("_in_" + name) )
        })
    }

    private void visitProcessOutputs(Statement outputs) {
        asDirectives(outputs).forEach((call) -> {
            fixLazyGString(call)

            final name = call.getMethodAsString()
            varToConstX(call.getArguments(), "tuple".equals(name), "each".equals(name))
            call.setMethod( constX("_out_" + name) )
            visitProcessOutputEmitAndTopic(call)
        })
    }

    private static final List<String> EMIT_AND_TOPIC = List.of("emit", "topic")

    private void visitProcessOutputEmitAndTopic(MethodCallExpression output) {
        final namedArgs = asNamedArgs(output)
        for( int i = 0; i < namedArgs.size(); i++ ) {
            final entry = namedArgs.get(i)
            final key = asConstX(entry.getKeyExpression())
            final value = asVarX(entry.getValueExpression())
            if( value != null && key != null && EMIT_AND_TOPIC.contains(key.getText()) ) {
                namedArgs.set(i, entryX(key, constX(value.getText())))
            }
        }
    }

    private void fixLazyGString(Expression node) {
        new GStringToLazyVisitor(sourceUnit).visit(node)
    }

    private Expression varToConstX(Expression node, boolean withinTuple, boolean withinEach) {
        if( node instanceof TupleExpression ) {
            final arguments = node.getExpressions()
            for( int i = 0; i < arguments.size(); i++ )
                arguments.set(i, varToConstX(arguments.get(i), withinTuple, withinEach))
            return node
        }

        if( node instanceof VariableExpression ) {
            final name = node.getName()

            if( "stdin".equals(name) && withinTuple )
                return createX( TokenStdinCall.class )

            if ( "stdout".equals(name) && withinTuple )
                return createX( TokenStdoutCall.class )

            return createX( TokenVar.class, constX(name) )
        }

        if( node instanceof MethodCallExpression ) {
            final name = node.getMethodAsString()
            final arguments = node.getArguments()

            if( "env".equals(name) && withinTuple )
                return createX( TokenEnvCall.class, (TupleExpression) varToStrX(arguments) )

            if( "eval".equals(name) && withinTuple )
                return createX( TokenEvalCall.class, (TupleExpression) varToStrX(arguments) )

            if( "file".equals(name) && (withinTuple || withinEach) )
                return createX( TokenFileCall.class, (TupleExpression) varToConstX(arguments, withinTuple, withinEach) )

            if( "path".equals(name) && (withinTuple || withinEach) )
                return createX( TokenPathCall.class, (TupleExpression) varToConstX(arguments, withinTuple, withinEach) )

            if( "val".equals(name) && withinTuple )
                return createX( TokenValCall.class, (TupleExpression) varToStrX(arguments) )
        }

        if( node instanceof PropertyExpression ) {
            // before:
            //   val( x.foo )
            // after:
            //   val({ x.foo })
            return wrapExpressionInClosure(node)
        }

        return node
    }

    private Expression varToStrX(Expression node) {
        if( node instanceof TupleExpression ) {
            final arguments = node.getExpressions()
            for( int i = 0; i < arguments.size(); i++ )
                arguments.set(i, varToStrX(arguments.get(i)))
            return node
        }

        if( node instanceof VariableExpression ) {
            // before:
            //   val(x)
            // after:
            //   val(TokenVar('x'))
            final name = node.getName()
            return createX( TokenVar.class, constX(name) )
        }

        if( node instanceof PropertyExpression ) {
            // before:
            //   tuple val( x.foo )
            // after:
            //   tuple val({ x.foo })
            return wrapExpressionInClosure(node)
        }

        return node
    }

    protected ClosureExpression wrapExpressionInClosure(Expression node)  {
        return closureX(block(new VariableScope(), stmt(node)))
    }

    private Statement processWhen(Expression when) {
        if( when instanceof EmptyExpression )
            return EmptyStatement.INSTANCE
        return stmt(callThisX("when", createX(
            TaskClosure.class,
            args(
                wrapExpressionInClosure(when),
                constX(getSourceText(when))
            )
        )))
    }

    private Statement processStub(Statement stub) {
        if( stub instanceof EmptyStatement )
            return EmptyStatement.INSTANCE
        return stmt(callThisX("stub", createX(
            TaskClosure.class,
            args(
                closureX(stub),
                constX(getSourceText(stub))
            )
        )))
    }

    @Override
    public void visitFunction(FunctionNode node) {
        if( RESERVED_NAMES.contains(node.getName()) ) {
            syntaxError(node, "`${node.getName()}` is not allowed as a function name because it is reserved for internal use")
            return
        }
        moduleNode.getScriptClassDummy().addMethod(node)
    }

    @Override
    public void visitOutput(OutputNode node) {
        visitOutputTargets(node.body)

        final closure = closureX(node.body)
        final result = stmt(callThisX("output", args(closure)))
        moduleNode.addStatement(result)
    }

    private void visitOutputTargets(Statement body) {
        for( final stmt : asBlockStatements(body) ) {
            final es = (ExpressionStatement)stmt
            final mce = (MethodCallExpression)es.getExpression()
            final name = mce.getMethod()
            final targetArgs = (ArgumentListExpression)mce.getArguments()
            final targetBody = (ClosureExpression)targetArgs[0]
            es.setExpression( callThisX('target', args(name, targetBody)) )
        }
    }

    // see: VariableScopeVisitor::visitExpressionStatement()
    @Override
    public void visitExpressionStatement(ExpressionStatement node) {
        final exp = node.getExpression()
        if( exp instanceof DeclarationExpression && exp.getNodeMetaData(IMPLICIT_DECLARATION) ) {
            final result = new AssignmentExpression(exp.getLeftExpression(), exp.getRightExpression())
            result.setSourcePosition(exp)
            node.setExpression(result)
        }
    }

    private String getSourceText(Statement node) {
        final builder = new StringBuilder()
        final colx = node.getColumnNumber()
        final colz = node.getLastColumnNumber()
        final first = node.getLineNumber()
        final last = node.getLastLineNumber()
        for( int i = first; i <= last; i++ ) {
            final line = sourceUnit.getSource().getLine(i, null)

            // prepend first-line indent
            if( i == first ) {
                int k = 0
                while( k < line.length() && line[k] == ' ' )
                    k++
                builder.append( line.substring(0, k) )
            }

            final begin = (i == first) ? colx - 1 : 0
            final end = (i == last) ? colz - 1 : line.length()
            builder.append( line.substring(begin, end) ).append('\n')
        }
        return builder.toString()
    }

    private String getSourceText(Expression node) {
        final stm = stmt(node)
        stm.setSourcePosition(node)
        return getSourceText(stm)
    }

    private String getSourceText(WorkflowNode node) {
        if( node.isEntry() && node.getLineNumber() == -1 )
            return getSourceText(node.main)

        final builder = new StringBuilder()
        final colx = node.getColumnNumber()
        final colz = node.getLastColumnNumber()
        final first = node.getLineNumber()
        final last = node.getLastLineNumber()
        for( int i = first; i <= last; i++ ) {
            def line = sourceUnit.getSource().getLine(i, null)
            if( i == last ) {
                line = line.substring(0, colz-1).replaceFirst(/}.*$/, '')
                if( !line.trim() )
                    continue
            }
            if( i == first ) {
                line = line.substring(colx-1).replaceFirst(/^.*\{/, '').trim()
                if( !line )
                    continue
            }
            builder.append(line).append('\n')
        }
        return builder.toString()
    }

    private void syntaxError(ASTNode node, String message) {
        sourceUnit.addError(new SyntaxException(message, node))
    }

    private static final String IMPLICIT_DECLARATION = "_IMPLICIT_DECLARATION"

}
