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

import java.lang.reflect.Modifier;
import java.util.ArrayList;
import java.util.Arrays;

import nextflow.script.ast.ASTNodeMarker;
import nextflow.script.ast.AgentNode;
import nextflow.script.ast.AssignmentExpression;
import nextflow.script.ast.RecordNode;
import nextflow.script.ast.ScriptNode;
import nextflow.script.ast.TupleParameter;
import org.codehaus.groovy.ast.ClassHelper;
import org.codehaus.groovy.ast.ClassNode;
import org.codehaus.groovy.ast.CodeVisitorSupport;
import org.codehaus.groovy.ast.FieldNode;
import org.codehaus.groovy.ast.Parameter;
import org.codehaus.groovy.ast.VariableScope;
import org.codehaus.groovy.ast.expr.ArgumentListExpression;
import org.codehaus.groovy.ast.expr.Expression;
import org.codehaus.groovy.ast.expr.MethodCallExpression;
import org.codehaus.groovy.ast.expr.VariableExpression;
import org.codehaus.groovy.ast.stmt.BlockStatement;
import org.codehaus.groovy.ast.stmt.ExpressionStatement;
import org.codehaus.groovy.ast.stmt.Statement;
import org.codehaus.groovy.control.SourceUnit;

import static nextflow.script.ast.ASTUtils.*;
import static nextflow.script.types.TypeCheckingUtils.getType;
import static org.codehaus.groovy.ast.tools.GeneralUtils.*;

/**
 * Lowers an {@link AgentNode} to a runtime {@code agent('name', { ... })} call.
 * The generated closure carries the directives, typed inputs/outputs and a
 * {@code PromptDef}, mirroring how {@link ProcessToGroovyVisitorV2} lowers a
 * process body.
 */
public class AgentToGroovyVisitor {

    private SourceUnit sourceUnit;

    private ScriptToGroovyHelper sgh;

    public AgentToGroovyVisitor(SourceUnit sourceUnit) {
        this.sourceUnit = sourceUnit;
        this.sgh = new ScriptToGroovyHelper(sourceUnit);
    }

    public Statement transform(AgentNode node) {
        // an agent's typed I/O IS a process's typed I/O, so the implicit stagers and the
        // output unstagers are inferred by the very same compiler units the process uses
        var stagers = new BlockStatement();
        for( var input : asFlatParams(node.inputs) )
            ImplicitStagers.visitInputType(input, varX(input.getName()), stagers);

        var unstagers = new BlockStatement();
        // one visitor for the whole agent, so the `$path<n>` keys are unique across outputs;
        // filesOnly because an agent has no task script to read `env`/`eval` back from
        var unstageVisitor = new ProcessToGroovyVisitorV2.ProcessUnstageVisitor(unstagers, true);
        var stdoutVisitor = new AgentStdoutVisitor();
        for( var stmt : asBlockStatements(node.outputs) ) {
            unstageVisitor.visit(stmt);
            // infer the output type before `stdout()` is rewritten
            getType(((ExpressionStatement) stmt).getExpression());
            stmt.visit(stdoutVisitor);
        }

        var statements = new ArrayList<Statement>();
        statements.add(node.directives);
        // the stagers/unstagers MUST precede the prompt: BaseScript.agent takes the closure's
        // RETURN value as the PromptDef, so the prompt statement has to stay last
        statements.add(stagers);
        statements.add(unstagers);
        statements.add(agentInputs(node.inputs));
        statements.add(agentOutputs(node.getName(), node.outputs));
        statements.add(agentPrompt(node.prompt));
        var body = closureX(block(new VariableScope(), statements));
        return stmt(callThisX("agent", args(constX(node.getName()), body)));
    }

    private Statement agentInputs(Parameter[] inputs) {
        var statements = Arrays.stream(inputs)
            .map((input) -> {
                var type = input.getType();
                if( input instanceof TupleParameter tp ) {
                    var components = Arrays.stream(tp.components)
                        .map(ProcessToGroovyVisitorV2::processInputCtor)
                        .toList();
                    return (Statement) stmt(callThisX("_input_", args(listX(components), classX(type))));
                }
                // a `Path?` declaration must mean optional here exactly as it does for a process,
                // otherwise a null value is rejected by TaskProcessor telling the user to append
                // the `?` they already appended
                var optional = type.getNodeMetaData(ASTNodeMarker.NULLABLE) != null;
                return (Statement) stmt(callThisX("_input_", args(constX(input.getName()), classX(type), constX(optional))));
            })
            .toList();
        return block(null, statements);
    }

    private Statement agentOutputs(String agentName, Statement outputs) {
        var statements = asBlockStatements(outputs).stream()
            .map(s -> ((ExpressionStatement) s).getExpression())
            .map(output -> outputDeclaration(agentName, output))
            .toList();
        return block(null, statements);
    }

    /**
     * One `output:` entry lowered to its `_output_` call. A bare variable is answered by the
     * model; an explicit right-hand side or a bare expression IS the output's value -- the
     * process rule verbatim -- so the model is neither asked for it nor allowed to bind it.
     */
    private Statement outputDeclaration(String agentName, Expression output) {
        if( output instanceof VariableExpression ve )
            return stmt(callThisX("_output_", outputType(agentName, ve.getName(), ve.getType())));
        if( output instanceof AssignmentExpression ae && ae.getLeftExpression() instanceof VariableExpression ve ) {
            var arguments = outputType(agentName, ve.getName(), ve.getType());
            arguments.addExpression(closureX(stmt(ae.getRightExpression())));
            return stmt(callThisX("_output_", arguments));
        }
        return stmt(callThisX("_output_", args(constX("$out"), classX(ProcessToGroovyVisitorV2.outputType(output)), closureX(stmt(output)))));
    }

    // A parameterized output type (e.g. `List<Integer>`) is preserved
    // from type erasure by storing it in a field of a hidden class, so
    // that the element type is available at runtime (see also
    // StripTypesVisitor).
    private ArgumentListExpression outputType(String agentName, String name, ClassNode type) {
        if( type.getGenericsTypes() == null )
            return args(constX(name), classX(type));
        var moduleNode = (ScriptNode) sourceUnit.getAST();
        var outputType = new RecordNode(ScriptToGroovyHelper.packageName(moduleNode) + "." + "__Output_" + agentName);
        outputType.addField(new FieldNode(name, Modifier.PUBLIC, type, outputType, null));
        moduleNode.addClass(outputType);
        return args(constX(name), classX(outputType), constX(name));
    }

    /**
     * Rewrite `stdout()` to `AgentOutputPlan.answer(stdout())`, since the task
     * stdout holds the terminal result frame rather than the answer itself.
     */
    private static class AgentStdoutVisitor extends CodeVisitorSupport {

        @Override
        public void visitMethodCallExpression(MethodCallExpression node) {
            if( node.isImplicitThis() && "stdout".equals(node.getMethodAsString()) && asMethodCallArguments(node).isEmpty() ) {
                node.setObjectExpression(classX(ClassHelper.make("nextflow.agent.AgentOutputPlan")));
                node.setMethod(constX("answer"));
                node.setArguments(args(callThisX("stdout")));
                node.setImplicitThis(false);
                return;
            }
            super.visitMethodCallExpression(node);
        }
    }

    private Statement agentPrompt(Statement prompt) {
        // the prompt is a block, exactly like a process body: the closure's value is its
        // last expression, so helper statements may precede the prompt text
        return stmt(createX(
            "nextflow.script.PromptDef",
            args(
                closureX(prompt),
                constX(sgh.getSourceText(prompt)),
                // capture the prompt closure's free-variable refs (params.*, task.ext.*)
                // so prompt-globals fold into the resume cache key (design §7.2/D3);
                // reuses the exact collector that populates process-body BodyDef.valRefs
                sgh.getVariableRefs(prompt)
            )
        ));
    }
}
