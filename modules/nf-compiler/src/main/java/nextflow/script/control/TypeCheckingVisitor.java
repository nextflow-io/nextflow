/*
 * Copyright 2024-2025, Seqera Labs
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

import java.util.Optional;

import nextflow.script.ast.ASTNodeMarker;
import nextflow.script.ast.AssignmentExpression;
import nextflow.script.ast.FeatureFlagNode;
import nextflow.script.ast.ProcessNode;
import nextflow.script.ast.ScriptNode;
import nextflow.script.ast.ScriptVisitorSupport;
import nextflow.script.ast.WorkflowNode;
import nextflow.script.types.Types;
import org.codehaus.groovy.ast.ASTNode;
import org.codehaus.groovy.ast.MethodNode;
import org.codehaus.groovy.ast.Variable;
import org.codehaus.groovy.ast.expr.Expression;
import org.codehaus.groovy.ast.expr.MethodCallExpression;
import org.codehaus.groovy.ast.expr.PropertyExpression;
import org.codehaus.groovy.ast.expr.VariableExpression;
import org.codehaus.groovy.ast.stmt.ExpressionStatement;
import org.codehaus.groovy.control.SourceUnit;
import org.codehaus.groovy.control.messages.SyntaxErrorMessage;
import org.codehaus.groovy.syntax.SyntaxException;

import static nextflow.script.ast.ASTHelpers.*;

/**
 * Validate types where possible.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public class TypeCheckingVisitor extends ScriptVisitorSupport {

    private SourceUnit sourceUnit;

    private boolean experimental;

    public TypeCheckingVisitor(SourceUnit sourceUnit, boolean experimental) {
        this.sourceUnit = sourceUnit;
        this.experimental = experimental;
    }

    @Override
    protected SourceUnit getSourceUnit() {
        return sourceUnit;
    }

    public void visit() {
        var moduleNode = sourceUnit.getAST();
        if( moduleNode instanceof ScriptNode sn )
            visit(sn);
    }

    // script declarations

    @Override
    public void visitFeatureFlag(FeatureFlagNode node) {
        if( !experimental )
            return;
        var fn = node.target;
        if( fn == null )
            return;
        var expectedType = fn.getType();
        var actualType = node.value.getType();
        if( !Types.isAssignableFrom(expectedType, actualType) )
            addError("Type mismatch for feature flag '" + node.name + "' -- expected a " + Types.getName(expectedType) + " but received a " + Types.getName(actualType), node);
    }

    // expressions

    @Override
    public void visitMethodCallExpression(MethodCallExpression node) {
        var defNode = (MethodNode) node.getNodeMetaData(ASTNodeMarker.METHOD_TARGET);
        if( defNode instanceof ProcessNode || defNode instanceof WorkflowNode )
            checkMethodCallArguments(node, defNode);
        super.visitMethodCallExpression(node);
    }

    private void checkMethodCallArguments(MethodCallExpression node, MethodNode defNode) {
        var argsCount = asMethodCallArguments(node).size();
        var paramsCount = numberOfParameters(defNode);
        if( argsCount != paramsCount )
            addError(String.format("Incorrect number of call arguments, expected %d but received %d", paramsCount, argsCount), node);
    }

    private static int numberOfParameters(MethodNode node) {
        if( node instanceof ProcessNode pn ) {
            return (int) asBlockStatements(pn.inputs).size();
        }
        if( node instanceof WorkflowNode wn ) {
            return (int) asBlockStatements(wn.takes).size();
        }
        return node.getParameters().length;
    }

    @Override
    public void visitPropertyExpression(PropertyExpression node) {
        var mn = asMethodOutput(node);
        if( mn instanceof ProcessNode pn )
            checkProcessOut(node, pn);
        else if( mn instanceof WorkflowNode wn )
            checkWorkflowOut(node, wn);
        else
            super.visitPropertyExpression(node);
    }

    private MethodNode asMethodOutput(PropertyExpression node) {
        if( node.getObjectExpression() instanceof PropertyExpression pe ) {
            if( pe.getObjectExpression() instanceof VariableExpression ve && "out".equals(pe.getPropertyAsString()) )
                return asMethodVariable(ve.getAccessedVariable());
        }
        return null;
    }

    private void checkProcessOut(PropertyExpression node, ProcessNode process) {
        var property = node.getPropertyAsString();
        var result = asDirectives(process.outputs)
            .filter(call -> property.equals(getProcessEmitName(call)))
            .findFirst();

        if( !result.isPresent() )
            addError("Unrecognized output `" + property + "` for process `" + process.getName() + "`", node);
    }

    private String getProcessEmitName(MethodCallExpression output) {
        return Optional.of(output)
            .flatMap(call -> Optional.ofNullable(asNamedArgs(call)))
            .flatMap(namedArgs ->
                namedArgs.stream()
                    .filter(entry -> "emit".equals(entry.getKeyExpression().getText()))
                    .findFirst()
            )
            .flatMap(entry -> Optional.ofNullable(
                entry.getValueExpression() instanceof VariableExpression ve ? ve.getName() : null
            ))
            .orElse(null);
    }

    private void checkWorkflowOut(PropertyExpression node, WorkflowNode workflow) {
        var property = node.getPropertyAsString();
        var result = asBlockStatements(workflow.emits).stream()
            .map(stmt -> ((ExpressionStatement) stmt).getExpression())
            .filter(emit -> property.equals(getWorkflowEmitName(emit)))
            .findFirst();

        if( !result.isPresent() )
            addError("Unrecognized output `" + property + "` for workflow `" + workflow.getName() + "`", node);
    }

    private String getWorkflowEmitName(Expression emit) {
        if( emit instanceof VariableExpression ve ) {
            return ve.getName();
        }
        else if( emit instanceof AssignmentExpression ae ) {
            var left = (VariableExpression)ae.getLeftExpression();
            return left.getName();
        }
        return null;
    }

    @Override
    public void addError(String message, ASTNode node) {
        var cause = new TypeError(message, node);
        var errorMessage = new SyntaxErrorMessage(cause, sourceUnit);
        sourceUnit.getErrorCollector().addErrorAndContinue(errorMessage);
    }

    private class TypeError extends SyntaxException implements PhaseAware {

        public TypeError(String message, ASTNode node) {
            super(message, node);
        }

        @Override
        public int getPhase() {
            return Phases.TYPE_CHECKING;
        }
    }
}
