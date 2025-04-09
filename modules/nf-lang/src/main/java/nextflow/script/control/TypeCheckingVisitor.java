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

import nextflow.script.ast.ASTNodeMarker;
import nextflow.script.ast.AssignmentExpression;
import nextflow.script.ast.FeatureFlagNode;
import nextflow.script.ast.ProcessNode;
import nextflow.script.ast.ScriptNode;
import nextflow.script.ast.ScriptVisitorSupport;
import nextflow.script.ast.WorkflowNode;
import nextflow.script.types.TypeChecker;
import nextflow.script.types.Types;
import org.codehaus.groovy.ast.ASTNode;
import org.codehaus.groovy.ast.MethodNode;
import org.codehaus.groovy.ast.Variable;
import org.codehaus.groovy.ast.expr.BinaryExpression;
import org.codehaus.groovy.ast.expr.DeclarationExpression;
import org.codehaus.groovy.ast.expr.Expression;
import org.codehaus.groovy.ast.expr.MethodCallExpression;
import org.codehaus.groovy.ast.expr.PropertyExpression;
import org.codehaus.groovy.ast.expr.VariableExpression;
import org.codehaus.groovy.ast.stmt.ExpressionStatement;
import org.codehaus.groovy.control.SourceUnit;
import org.codehaus.groovy.control.messages.SyntaxErrorMessage;
import org.codehaus.groovy.syntax.SyntaxException;

import static nextflow.script.ast.ASTUtils.*;

/**
 * Resolve and validate the types of expressions.
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

    // statements

    @Override
    public void visitExpressionStatement(ExpressionStatement node) {
        var exp = node.getExpression();
        if( exp instanceof AssignmentExpression ae && ae.getNodeMetaData(ASTNodeMarker.IMPLICIT_DECLARATION) != null ) {
            applyDeclaration(ae);
            return;
        }
        super.visitExpressionStatement(node);
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
    public void visitDeclarationExpression(DeclarationExpression node) {
        applyDeclaration(node);
    }

    private void applyDeclaration(BinaryExpression node) {
        var target = node.getLeftExpression();
        var source = node.getRightExpression();
        visit(target);
        visit(source);
        var sourceType = TypeChecker.getType(source);
        target.putNodeMetaData(ASTNodeMarker.INFERRED_TYPE, sourceType);
    }

    @Override
    public void visitPropertyExpression(PropertyExpression node) {
        super.visitPropertyExpression(node);

        if( TypeChecker.getType(node) == null ) {
            var mn = asMethodNamedOutput(node);
            var property = node.getPropertyAsString();
            if( mn instanceof ProcessNode pn )
                addError("Unrecognized output `" + property + "` for process `" + pn.getName() + "`", node);
            else if( mn instanceof WorkflowNode wn )
                addError("Unrecognized output `" + property + "` for workflow `" + wn.getName() + "`", node);
        }
    }

    private static MethodNode asMethodNamedOutput(PropertyExpression node) {
        if( node.getObjectExpression() instanceof PropertyExpression pe )
            return asMethodOutput(pe);
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
