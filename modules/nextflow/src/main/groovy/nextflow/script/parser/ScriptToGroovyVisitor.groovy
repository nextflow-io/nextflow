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

import nextflow.script.BaseScript
import nextflow.script.BodyDef
import nextflow.script.IncludeDef
import nextflow.script.TokenEnvCall
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
import org.codehaus.groovy.ast.expr.EmptyExpression
import org.codehaus.groovy.ast.expr.Expression
import org.codehaus.groovy.ast.expr.MethodCallExpression
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
 * Transform Nextflow AST nodes into pure Groovy.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
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
        final moduleArgs = node.modules.stream()
            .map((module) -> {
                final name = constX(module.name)
                return module.alias != null
                    ? (Expression) createX(IncludeDef.Module.class, name, constX(module.alias))
                    : (Expression) createX(IncludeDef.Module.class, name)
            })
            .collect(Collectors.toList())

        final include = callThisX("include", args(createX(IncludeDef.class, args((Expression) listX(moduleArgs), (Expression) node.source))))
        final result = stmt(callX(include, "load0", new ArgumentListExpression()))
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
        visitWorkflowEmits(node.emits, node.main)
        visitWorkflowPublishers(node.publishers, node.main)

        final bodyDef = stmt(createX(
            BodyDef.class,
            args(
                closureX(node.main),
                constX(""), // TODO: source code formatting
                constX("workflow")
            )
        ))
        final closure = closureX(block(new VariableScope(), List.of(
            node.takes, // TODO: cannot be null
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
                final target = varX("$out")
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
        visitProcessInputs(node.inputs)
        visitProcessOutputs(node.outputs)

        final when = processWhen(node.when)
        final bodyDef = stmt(createX(
            BodyDef.class,
            args(
                closureX(node.exec),
                constX(""), // TODO: source code formatting
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

    private void visitProcessInputs(Statement inputs) {
        asDirectives(inputs).forEach((call) -> {
            final name = call.getMethodAsString()
            varToConstX(call.getArguments(), "tuple".equals(name), "each".equals(name))
            call.setMethod( constX("_in_" + name) )
        })
    }

    private void visitProcessOutputs(Statement outputs) {
        asDirectives(outputs).forEach((call) -> {
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

    private Expression varToConstX(Expression node, boolean withinTuple, boolean withinEach) {
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

            if( "file".equals(name) && (withinTuple || withinEach) )
                return createX( TokenFileCall.class, (TupleExpression) varToConstX(arguments, withinTuple, withinEach) )

            if( "path".equals(name) && (withinTuple || withinEach) )
                return createX( TokenPathCall.class, (TupleExpression) varToConstX(arguments, withinTuple, withinEach) )

            if( "val".equals(name) && withinTuple )
                return createX( TokenValCall.class, (TupleExpression) varToStrX(arguments) )
        }

        if( node instanceof TupleExpression ) {
            final arguments = node.getExpressions()
            for( int i = 0; i < arguments.size(); i++ )
                arguments.set(i, varToConstX(arguments.get(i), withinTuple, withinEach))
            return node
        }

        return node
    }

    private Expression varToStrX(Expression node) {
        if( node instanceof VariableExpression ) {
            final name = node.getName()
            return createX( TokenVar.class, constX(name) )
        }

        if( node instanceof TupleExpression ) {
            final arguments = node.getExpressions()
            for( int i = 0; i < arguments.size(); i++ )
                arguments.set(i, varToStrX(arguments.get(i)))
            return node
        }

        return node
    }

    private Statement processWhen(Expression when) {
        if( when instanceof EmptyExpression )
            return EmptyStatement.INSTANCE
        return stmt(when)
    }

    private Statement processStub(Statement stub) {
        if( stub instanceof EmptyStatement )
            return EmptyStatement.INSTANCE
        return stmt(callThisX("stub", closureX(stub)))
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

    private void syntaxError(ASTNode node, String message) {
        sourceUnit.addError(new SyntaxException(message, node))
    }

}
