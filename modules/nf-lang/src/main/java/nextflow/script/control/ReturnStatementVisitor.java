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

import java.util.ArrayList;
import java.util.List;

import nextflow.script.ast.AssignmentExpression;
import nextflow.script.dsl.Types;
import org.codehaus.groovy.ast.ASTNode;
import org.codehaus.groovy.ast.ClassHelper;
import org.codehaus.groovy.ast.ClassNode;
import org.codehaus.groovy.ast.ClassCodeVisitorSupport;
import org.codehaus.groovy.ast.expr.ConstantExpression;
import org.codehaus.groovy.ast.expr.DeclarationExpression;
import org.codehaus.groovy.ast.expr.Expression;
import org.codehaus.groovy.ast.stmt.BlockStatement;
import org.codehaus.groovy.ast.stmt.CatchStatement;
import org.codehaus.groovy.ast.stmt.ExpressionStatement;
import org.codehaus.groovy.ast.stmt.IfStatement;
import org.codehaus.groovy.ast.stmt.ReturnStatement;
import org.codehaus.groovy.ast.stmt.Statement;
import org.codehaus.groovy.ast.stmt.ThrowStatement;
import org.codehaus.groovy.ast.stmt.TryCatchStatement;
import org.codehaus.groovy.control.ErrorCollector;
import org.codehaus.groovy.control.SourceUnit;
import org.codehaus.groovy.control.messages.SyntaxErrorMessage;

import static nextflow.script.ast.ASTUtils.*;
import static nextflow.script.types.TypeCheckingUtils.*;

/**
 * Infer the return type of a code block based on implicit
 * and explicit return statements.
 *
 * @see org.codehaus.groovy.classgen.ReturnAdder
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public class ReturnStatementVisitor extends ClassCodeVisitorSupport {

    private SourceUnit sourceUnit;

    private ErrorCollector errorCollector;

    private ClassNode returnType;

    private ClassNode inferredReturnType;

    private boolean coerce;

    private List<ASTNode> missingReturns = new ArrayList<>();

    public ReturnStatementVisitor(SourceUnit sourceUnit, ErrorCollector errorCollector) {
        this.sourceUnit = sourceUnit;
        this.errorCollector = errorCollector;
    }

    @Override
    protected SourceUnit getSourceUnit() {
        return sourceUnit;
    }

    public void visit(ASTNode owner, ClassNode returnType, Statement code) {
        visit(owner, returnType, code, false);
    }

    /**
     * @param owner      function or closure that contains the code
     * @param returnType
     * @param code
     * @param coerce     accept any return value that can be coerced to the return type
     */
    public void visit(ASTNode owner, ClassNode returnType, Statement code, boolean coerce) {
        this.returnType = returnType;
        this.coerce = coerce;
        visit(addReturnsIfNeeded(code, owner));
        if( returnsValue() ) {
            for( var node : missingReturns )
                addError("Missing return statement", node);
        }
        this.returnType = null;
    }

    /**
     * Convert trailing expression statements into return statements,
     * and record any code path that ends without a return value.
     *
     * @param node
     * @param parent  node to report if the statement has no source position
     */
    private Statement addReturnsIfNeeded(Statement node, ASTNode parent) {
        if( node instanceof BlockStatement block && !block.isEmpty() ) {
            var statements = new ArrayList<>(block.getStatements());
            int lastIndex = statements.size() - 1;
            var last = addReturnsIfNeeded(statements.get(lastIndex), block);
            statements.set(lastIndex, last);
            return withSourcePosition(new BlockStatement(statements, block.getVariableScope()), block);
        }

        if( node instanceof ExpressionStatement es && !isAssignment(es.getExpression()) ) {
            return withSourcePosition(new ReturnStatement(es.getExpression()), es);
        }

        if( node instanceof IfStatement ies ) {
            return withSourcePosition(new IfStatement(
                ies.getBooleanExpression(),
                addReturnsIfNeeded(ies.getIfBlock(), ies),
                addReturnsIfNeeded(ies.getElseBlock(), ies) ), ies);
        }

        if( node instanceof TryCatchStatement tcs ) {
            var result = new TryCatchStatement(addReturnsIfNeeded(tcs.getTryStatement(), tcs), tcs.getFinallyStatement());
            for( var cs : tcs.getCatchStatements() )
                result.addCatch(withSourcePosition(new CatchStatement(cs.getVariable(), addReturnsIfNeeded(cs.getCode(), cs)), cs));
            return withSourcePosition(result, tcs);
        }

        if( !(node instanceof ReturnStatement) && !(node instanceof ThrowStatement) )
            missingReturns.add(node.getLineNumber() != -1 && !(node instanceof BlockStatement) ? node : parent);

        return node;
    }

    private static boolean isAssignment(Expression node) {
        return node instanceof DeclarationExpression || node instanceof AssignmentExpression;
    }

    private boolean returnsValue() {
        if( !ClassHelper.isDynamicTyped(returnType) )
            return !ClassHelper.VOID_TYPE.equals(returnType);
        return inferredReturnType != null && !ClassHelper.VOID_TYPE.equals(inferredReturnType);
    }

    private static <T extends Statement> T withSourcePosition(T node, Statement source) {
        node.setSourcePosition(source);
        return node;
    }

    @Override
    public void visitReturnStatement(ReturnStatement node) {
        // a bare `return` yields no value, unlike `return null`
        var expression = node.getExpression();
        var sourceType = expression == ConstantExpression.EMPTY_EXPRESSION
            ? ClassHelper.VOID_TYPE
            : getType(expression);
        if( coerce ) {
            // a void expression has no value to coerce
            if( ClassHelper.VOID_TYPE.equals(sourceType) )
                addError(String.format("Return value with type void does not match the declared return type (%s)", Types.getName(returnType)), node);
            return;
        }
        if( inferredReturnType != null && !ClassHelper.isDynamicTyped(returnType) ) {
            if( !Types.isAssignableFrom(inferredReturnType, sourceType) )
                addError(String.format("Return value with type %s does not match previous return type (%s)", Types.getName(sourceType), Types.getName(inferredReturnType)), node);
        }
        else if( Types.isAssignableFrom(returnType, sourceType) ) {
            inferredReturnType = sourceType;
        }
        else {
            addError(String.format("Return value with type %s does not match the declared return type (%s)", Types.getName(sourceType), Types.getName(returnType)), node);
        }
    }

    public ClassNode getInferredReturnType() {
        return inferredReturnType;
    }

    @Override
    public void addError(String message, ASTNode node) {
        var cause = new TypeError(message, node);
        var errorMessage = new SyntaxErrorMessage(cause, sourceUnit);
        errorCollector.addErrorAndContinue(errorMessage);
    }
}
