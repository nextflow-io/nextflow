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
import java.util.Arrays;
import java.util.List;

import nextflow.script.ast.ASTNodeMarker;
import nextflow.script.ast.AssignmentExpression;
import nextflow.script.ast.RecordNode;
import nextflow.script.ast.ScriptNode;
import nextflow.script.ast.WorkflowNode;
import org.codehaus.groovy.ast.FieldNode;
import org.codehaus.groovy.ast.Parameter;
import org.codehaus.groovy.ast.VariableScope;
import org.codehaus.groovy.ast.expr.BinaryExpression;
import org.codehaus.groovy.ast.expr.VariableExpression;
import org.codehaus.groovy.ast.stmt.BlockStatement;
import org.codehaus.groovy.ast.stmt.ExpressionStatement;
import org.codehaus.groovy.ast.stmt.Statement;
import org.codehaus.groovy.control.SourceUnit;

import static nextflow.script.ast.ASTUtils.*;
import static org.codehaus.groovy.ast.tools.GeneralUtils.*;

/**
 * Transform a workflow definition into a Groovy AST.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public class WorkflowToGroovyVisitor {

    private SourceUnit sourceUnit;

    private ScriptNode moduleNode;

    public WorkflowToGroovyVisitor(SourceUnit sourceUnit) {
        this.sourceUnit = sourceUnit;
        this.moduleNode = (ScriptNode) sourceUnit.getAST();
    }

    public Statement transform(WorkflowNode node) {
        var main = node.main instanceof BlockStatement block ? block : new BlockStatement();
        visitWorkflowEmits(node.emits, main);
        visitWorkflowPublishers(node.publishers, main);
        visitWorkflowHandler(node.onComplete, "setOnComplete", main);
        visitWorkflowHandler(node.onError, "setOnError", main);

        var bodyDef = stmt(createX(
            "nextflow.script.BodyDef",
            args(
                closureX(null, main),
                constX(null),
                constX("workflow")
            )
        ));
        var closure = closureX(null, block(new VariableScope(), List.of(
            workflowTakes(node.getParameters(), node.isEntry() ? null : node.getName()),
            node.emits,
            bodyDef
        )));
        var arguments = node.isEntry()
            ? args(closure)
            : args(constX(node.getName()), closure);
        return stmt(callThisX("workflow", arguments));
    }

    // The declared type of each workflow take is preserved from type
    // erasure by storing it in a field of a hidden class, so that e.g.
    // the element type of a `Channel<Sample>` input is available at
    // runtime (see also StripTypesVisitor).
    private Statement workflowTakes(Parameter[] takes, String workflowName) {
        if( takes.length == 0 )
            return block(null, List.of());

        var takesType = new RecordNode(ScriptToGroovyHelper.packageName(moduleNode) + "." + "__Takes"
            + (workflowName != null ? "_" + workflowName : ""));
        for( var take : takes ) {
            takesType.addField(new FieldNode(
                take.getName(),
                Modifier.PUBLIC,
                take.getType(),
                takesType,
                null
            ));
        }
        moduleNode.addClass(takesType);

        var statements = Arrays.stream(takes)
            .map((take) -> {
                var optional = take.getType().getNodeMetaData(ASTNodeMarker.NULLABLE) != null;
                return stmt(callThisX("_take_", args(constX(take.getName()), classX(takesType), constX(optional))));
            })
            .toList();
        return block(null, statements);
    }

    private void visitWorkflowEmits(Statement emits, BlockStatement main) {
        for( var stmt : asBlockStatements(emits) ) {
            var es = (ExpressionStatement)stmt;
            var emit = es.getExpression();
            if( emit instanceof VariableExpression ve ) {
                es.setExpression(callThisX("_emit_", args(constX(ve.getName()))));
            }
            else if( emit instanceof AssignmentExpression ae ) {
                var target = (VariableExpression)ae.getLeftExpression();
                main.addStatement(assignS(target, emit));
                es.setExpression(callThisX("_emit_", args(constX(target.getName()))));
                main.addStatement(es);
            }
            else {
                var target = varX("$out");
                main.addStatement(assignS(target, emit));
                es.setExpression(callThisX("_emit_", args(constX(target.getName()))));
                main.addStatement(es);
            }
        }
    }

    private void visitWorkflowPublishers(Statement publishers, BlockStatement main) {
        for( var stmt : asBlockStatements(publishers) ) {
            var es = (ExpressionStatement)stmt;
            var publish = (BinaryExpression)es.getExpression();
            var target = asVarX(publish.getLeftExpression());
            es.setExpression(callThisX("_publish_", args(constX(target.getName()), publish.getRightExpression())));
            main.addStatement(es);
        }
    }

    private void visitWorkflowHandler(Statement code, String name, BlockStatement main) {
        if( code instanceof BlockStatement block )
            main.addStatement(stmt(callX(varX("workflow"), name, args(closureX(null, block)))));
    }

}
