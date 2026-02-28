/*
 * Copyright 2013-2026, Seqera Labs
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
package nextflow.script.control;

import java.util.List;

import nextflow.script.ast.ProcessNodeV1;
import org.codehaus.groovy.ast.VariableScope;
import org.codehaus.groovy.ast.expr.ClosureExpression;
import org.codehaus.groovy.ast.expr.EmptyExpression;
import org.codehaus.groovy.ast.expr.Expression;
import org.codehaus.groovy.ast.expr.MethodCallExpression;
import org.codehaus.groovy.ast.expr.TupleExpression;
import org.codehaus.groovy.ast.expr.VariableExpression;
import org.codehaus.groovy.ast.stmt.EmptyStatement;
import org.codehaus.groovy.ast.stmt.Statement;
import org.codehaus.groovy.control.SourceUnit;

import static nextflow.script.ast.ASTUtils.*;
import static org.codehaus.groovy.ast.tools.GeneralUtils.*;

/**
 * Transform a legacy process AST node into Groovy AST.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public class ProcessToGroovyVisitorV1 {

    private SourceUnit sourceUnit;

    private ScriptToGroovyHelper sgh;

    public ProcessToGroovyVisitorV1(SourceUnit sourceUnit) {
        this.sourceUnit = sourceUnit;
        this.sgh = new ScriptToGroovyHelper(sourceUnit);
    }

    public Statement transform(ProcessNodeV1 node) {
        visitProcessDirectives(node.directives);
        visitProcessInputs(node.inputs);
        visitProcessOutputs(node.outputs);

        if( "script".equals(node.type) )
            node.exec.visit(new TaskCmdXformVisitor(sourceUnit));
        node.stub.visit(new TaskCmdXformVisitor(sourceUnit));

        var closure = closureX(null, block(new VariableScope(), List.of(
            node.directives,
            node.inputs,
            node.outputs,
            processWhen(node.when),
            processStub(node.stub),
            stmt(createX(
                "nextflow.script.BodyDef",
                args(
                    closureX(null, node.exec),
                    constX(sgh.getSourceText(node.exec)),
                    constX(node.type),
                    sgh.getVariableRefs(node.exec)
                )
            ))
        )));
        return stmt(callThisX("process", args(constX(node.getName()), closure)));
    }

    private void visitProcessDirectives(Statement directives) {
        asDirectives(directives).forEach((call) -> {
            var arguments = asMethodCallArguments(call);
            if( arguments.size() != 1 )
                return;
            var firstArg = arguments.get(0);
            if( firstArg instanceof ClosureExpression )
                return;
            arguments.set(0, sgh.transformToLazy(firstArg));
        });
    }

    private void visitProcessInputs(Statement inputs) {
        asDirectives(inputs).forEach((call) -> {
            var name = call.getMethodAsString();
            varToConstX(call.getArguments(), "tuple".equals(name), "each".equals(name));
            call.setMethod( constX("_in_" + name) );
        });
    }

    private void visitProcessOutputs(Statement outputs) {
        asDirectives(outputs).forEach((call) -> {
            var name = call.getMethodAsString();
            varToConstX(call.getArguments(), "tuple".equals(name), false);
            call.setMethod( constX("_out_" + name) );
            visitProcessOutputEmitAndTopic(call);
        });
    }

    private static final List<String> EMIT_AND_TOPIC = List.of("emit", "topic");

    private void visitProcessOutputEmitAndTopic(MethodCallExpression output) {
        var namedArgs = asNamedArgs(output);
        for( int i = 0; i < namedArgs.size(); i++ ) {
            var entry = namedArgs.get(i);
            var key = asConstX(entry.getKeyExpression());
            var value = asVarX(entry.getValueExpression());
            if( value != null && key != null && EMIT_AND_TOPIC.contains(key.getText()) ) {
                namedArgs.set(i, entryX(key, constX(value.getText())));
            }
        }
    }

    private Expression varToConstX(Expression node, boolean withinTuple, boolean withinEach) {
        if( node instanceof TupleExpression te ) {
            var arguments = te.getExpressions();
            for( int i = 0; i < arguments.size(); i++ )
                arguments.set(i, varToConstX(arguments.get(i), withinTuple, withinEach));
            return te;
        }

        if( node instanceof VariableExpression ve ) {
            var name = ve.getName();

            if( "stdin".equals(name) && withinTuple )
                return createX( "nextflow.script.TokenStdinCall" );

            if ( "stdout".equals(name) && withinTuple )
                return createX( "nextflow.script.TokenStdoutCall" );

            return createX( "nextflow.script.TokenVar", constX(name) );
        }

        if( node instanceof MethodCallExpression mce ) {
            var name = mce.getMethodAsString();
            var arguments = mce.getArguments();

            if( "env".equals(name) && withinTuple )
                return createX( "nextflow.script.TokenEnvCall", (TupleExpression) varToStrX(arguments) );

            if( "eval".equals(name) && withinTuple )
                return createX( "nextflow.script.TokenEvalCall", (TupleExpression) varToStrX(arguments) );

            if( "file".equals(name) && (withinTuple || withinEach) )
                return createX( "nextflow.script.TokenFileCall", (TupleExpression) varToConstX(arguments, withinTuple, withinEach) );

            if( "path".equals(name) && (withinTuple || withinEach) )
                return createX( "nextflow.script.TokenPathCall", (TupleExpression) varToConstX(arguments, withinTuple, withinEach) );

            if( "val".equals(name) && withinTuple )
                return createX( "nextflow.script.TokenValCall", (TupleExpression) varToStrX(arguments) );
        }

        return sgh.transformToLazy(node);
    }

    private Expression varToStrX(Expression node) {
        if( node instanceof TupleExpression te ) {
            var arguments = te.getExpressions();
            for( int i = 0; i < arguments.size(); i++ )
                arguments.set(i, varToStrX(arguments.get(i)));
            return te;
        }

        if( node instanceof VariableExpression ve ) {
            // before:
            //   val(x)
            // after:
            //   val(TokenVar('x'))
            var name = ve.getName();
            return createX( "nextflow.script.TokenVar", constX(name) );
        }

        return sgh.transformToLazy(node);
    }

    private Statement processWhen(Expression when) {
        if( when instanceof EmptyExpression )
            return EmptyStatement.INSTANCE;
        return stmt(callThisX("when", createX(
            "nextflow.script.TaskClosure",
            args(
                closureX(null, block(stmt(when))),
                constX(sgh.getSourceText(when))
            )
        )));
    }

    private Statement processStub(Statement stub) {
        if( stub instanceof EmptyStatement )
            return EmptyStatement.INSTANCE;
        return stmt(callThisX("stub", createX(
            "nextflow.script.TaskClosure",
            args(
                closureX(null, stub),
                constX(sgh.getSourceText(stub))
            )
        )));
    }

}
