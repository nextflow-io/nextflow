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
import nextflow.script.ast.ProcessNode;
import nextflow.script.ast.ScriptNode;
import nextflow.script.ast.ScriptVisitorSupport;
import nextflow.script.ast.WorkflowNode;
import org.codehaus.groovy.ast.ASTNode;
import org.codehaus.groovy.ast.MethodNode;
import org.codehaus.groovy.ast.expr.MethodCallExpression;
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

    public TypeCheckingVisitor(SourceUnit sourceUnit) {
        this.sourceUnit = sourceUnit;
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
        var paramsCount = defNode.getParameters().length;
        if( argsCount != paramsCount )
            addError(String.format("Incorrect number of call arguments, expected %d but received %d", paramsCount, argsCount), node);
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
