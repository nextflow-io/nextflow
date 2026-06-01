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

import nextflow.script.ast.AgentNode;
import nextflow.script.ast.AssignmentExpression;
import org.codehaus.groovy.ast.Parameter;
import org.codehaus.groovy.ast.VariableScope;
import org.codehaus.groovy.ast.expr.Expression;
import org.codehaus.groovy.ast.expr.VariableExpression;
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
public class AgentToGroovyVisitorV2 {

    private SourceUnit sourceUnit;

    private ScriptToGroovyHelper sgh;

    public AgentToGroovyVisitorV2(SourceUnit sourceUnit) {
        this.sourceUnit = sourceUnit;
        this.sgh = new ScriptToGroovyHelper(sourceUnit);
    }

    public Statement transform(AgentNode node) {
        var statements = new ArrayList<Statement>();
        statements.add(node.directives);
        statements.add(agentInputs(node.inputs));
        statements.add(agentOutputs(node.outputs));
        statements.add(agentPrompt(node.prompt));
        var body = closureX(block(new VariableScope(), statements));
        return stmt(callThisX("agent", args(constX(node.getName()), body)));
    }

    private Statement agentInputs(Parameter[] inputs) {
        var statements = Arrays.stream(inputs)
            .map((input) -> stmt(callThisX("_input_", args(constX(input.getName()), classX(input.getType())))))
            .map(s -> (Statement) s)
            .toList();
        return block(null, statements);
    }

    private Statement agentOutputs(Statement outputs) {
        var statements = new ArrayList<Statement>();
        for( var stmt : asBlockStatements(outputs) ) {
            var output = ((ExpressionStatement) stmt).getExpression();
            var target = targetOf(output);
            if( target != null )
                statements.add(stmt(callThisX("_output_", args(constX(target.getName()), classX(target.getType())))));
        }
        return block(null, statements);
    }

    private VariableExpression targetOf(Expression output) {
        if( output instanceof VariableExpression ve )
            return ve;
        if( output instanceof AssignmentExpression ae && ae.getLeftExpression() instanceof VariableExpression ve )
            return ve;
        return null;
    }

    private Statement agentPrompt(Statement prompt) {
        var expr = ((ExpressionStatement) prompt).getExpression();
        return stmt(createX(
            "nextflow.script.PromptDef",
            args(
                closureX(stmt(expr)),
                constX(sgh.getSourceText(prompt))
            )
        ));
    }
}
