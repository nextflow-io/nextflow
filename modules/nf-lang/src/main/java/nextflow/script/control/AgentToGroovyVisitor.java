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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Objects;

import nextflow.script.ast.ASTNodeMarker;
import nextflow.script.ast.AgentNode;
import nextflow.script.ast.AssignmentExpression;
import org.codehaus.groovy.ast.Parameter;
import org.codehaus.groovy.ast.VariableScope;
import org.codehaus.groovy.ast.expr.Expression;
import org.codehaus.groovy.ast.expr.VariableExpression;
import org.codehaus.groovy.ast.stmt.BlockStatement;
import org.codehaus.groovy.ast.stmt.ExpressionStatement;
import org.codehaus.groovy.ast.stmt.Statement;
import org.codehaus.groovy.control.SourceUnit;

import static nextflow.script.ast.ASTUtils.*;
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
        // deliberately the RAW inputs, not asFlatParams: a tuple input is rejected upstream
        // (ScriptResolveVisitor) precisely because it declares no context slot to stage from,
        // and this keeps the stager loop provably in step with the `_input_` loop below
        for( var input : node.inputs )
            ImplicitStagers.visitInputType(input, varX(input.getName()), stagers);

        var unstagers = new BlockStatement();
        // one visitor for the whole agent, so the `$path<n>` keys are unique across outputs;
        // filesOnly because an agent has no task script to read `env`/`eval` back from
        var unstageVisitor = new ProcessToGroovyVisitorV2.ProcessUnstageVisitor(unstagers, true);
        for( var stmt : asBlockStatements(node.outputs) )
            unstageVisitor.visit(stmt);

        var statements = new ArrayList<Statement>();
        statements.add(node.directives);
        // the stagers/unstagers MUST precede the prompt: BaseScript.agent takes the closure's
        // RETURN value as the PromptDef, so the prompt statement has to stay last
        statements.add(stagers);
        statements.add(unstagers);
        statements.add(agentInputs(node.inputs));
        statements.add(agentOutputs(node.outputs));
        statements.add(agentPrompt(node.prompt));
        var body = closureX(block(new VariableScope(), statements));
        return stmt(callThisX("agent", args(constX(node.getName()), body)));
    }

    private Statement agentInputs(Parameter[] inputs) {
        var statements = Arrays.stream(inputs)
            .map((input) -> {
                var type = input.getType();
                // a `Path?` declaration must mean optional here exactly as it does for a process,
                // otherwise a null value is rejected by TaskProcessor telling the user to append
                // the `?` they already appended
                var optional = type.getNodeMetaData(ASTNodeMarker.NULLABLE) != null;
                return (Statement) stmt(callThisX("_input_", args(constX(input.getName()), classX(type), constX(optional))));
            })
            .toList();
        return block(null, statements);
    }

    private Statement agentOutputs(Statement outputs) {
        var statements = asBlockStatements(outputs).stream()
            .map(s -> ((ExpressionStatement) s).getExpression())
            .map(AgentToGroovyVisitor::outputDeclaration)
            .filter(Objects::nonNull)
            .toList();
        return block(null, statements);
    }

    /**
     * One `output:` entry lowered to its `_output_` call, or null for a form that is not an output
     * declaration. A bare variable is answered by the model; an explicit right-hand side IS the
     * output's value -- the process rule verbatim -- so the model is neither asked for it nor
     * allowed to bind it.
     */
    private static Statement outputDeclaration(Expression output) {
        if( output instanceof VariableExpression ve )
            return stmt(callThisX("_output_", args(constX(ve.getName()), classX(ve.getType()))));
        if( output instanceof AssignmentExpression ae && ae.getLeftExpression() instanceof VariableExpression ve )
            return stmt(callThisX("_output_", args(constX(ve.getName()), classX(ve.getType()), closureX(stmt(ae.getRightExpression())))));
        return null;
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
