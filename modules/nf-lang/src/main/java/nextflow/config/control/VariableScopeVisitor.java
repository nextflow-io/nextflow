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
package nextflow.config.control;

import java.util.ArrayList;
import java.util.List;
import java.util.Stack;

import nextflow.config.ast.ConfigApplyBlockNode;
import nextflow.config.ast.ConfigApplyNode;
import nextflow.config.ast.ConfigAssignNode;
import nextflow.config.ast.ConfigBlockNode;
import nextflow.config.ast.ConfigIncludeNode;
import nextflow.config.ast.ConfigNode;
import nextflow.config.ast.ConfigVisitorSupport;
import nextflow.config.dsl.ConfigDsl;
import nextflow.config.schema.SchemaNode;
import nextflow.script.ast.ASTNodeMarker;
import nextflow.script.control.VariableScopeChecker;
import nextflow.script.dsl.ProcessDsl;
import nextflow.script.dsl.ScriptDsl;
import org.codehaus.groovy.ast.ClassNode;
import org.codehaus.groovy.ast.DynamicVariable;
import org.codehaus.groovy.ast.Variable;
import org.codehaus.groovy.ast.VariableScope;
import org.codehaus.groovy.ast.expr.BinaryExpression;
import org.codehaus.groovy.ast.expr.ClosureExpression;
import org.codehaus.groovy.ast.expr.ConstantExpression;
import org.codehaus.groovy.ast.expr.DeclarationExpression;
import org.codehaus.groovy.ast.expr.Expression;
import org.codehaus.groovy.ast.expr.MethodCallExpression;
import org.codehaus.groovy.ast.expr.TupleExpression;
import org.codehaus.groovy.ast.expr.VariableExpression;
import org.codehaus.groovy.ast.stmt.BlockStatement;
import org.codehaus.groovy.ast.stmt.CatchStatement;
import org.codehaus.groovy.ast.stmt.ExpressionStatement;
import org.codehaus.groovy.control.SourceUnit;
import org.codehaus.groovy.syntax.Types;

/**
 * Initialize the variable scopes for an AST.
 *
 * @see org.codehaus.groovy.classgen.VariableScopeVisitor
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class VariableScopeVisitor extends ConfigVisitorSupport {

    private SourceUnit sourceUnit;

    private VariableScopeChecker vsc;

    private Stack<String> configScopes = new Stack<>();

    public VariableScopeVisitor(SourceUnit sourceUnit) {
        this.sourceUnit = sourceUnit;
        this.vsc = new VariableScopeChecker(sourceUnit, new ClassNode(ConfigDsl.class));
    }

    @Override
    protected SourceUnit getSourceUnit() {
        return sourceUnit;
    }

    public void visit() {
        var moduleNode = sourceUnit.getAST();
        if( moduleNode instanceof ConfigNode cn ) {
            super.visit(cn);
            vsc.checkUnusedVariables();
        }
    }

    @Override
    public void visitConfigApplyBlock(ConfigApplyBlockNode node) {
        configScopes.add(node.name);
        var names = currentConfigScopes();
        var option = SchemaNode.ROOT.getDslOption(names);
        if( option != null ) {
            vsc.pushScope(option.dsl());
            super.visitConfigApplyBlock(node);
            vsc.popScope();
        }
        else {
            // invalid config apply block is handled by ScriptAstBuilder
            // addError("Unrecognized config block '" + node.name + "'", node);
        }
        configScopes.pop();
    }

    @Override
    public void visitConfigApply(ConfigApplyNode node) {
        checkMethodCall(node);
    }

    @Override
    public void visitConfigAssign(ConfigAssignNode node) {
        for( int i = 0; i < node.names.size() - 1; i++ )
            configScopes.add(node.names.get(i));

        var scopes = currentConfigScopes();
        var inProcess = !scopes.isEmpty() && "process".equals(scopes.get(0));
        var inClosure = node.value instanceof ClosureExpression;
        if( inClosure && !inProcess && !isWorkflowHandler(scopes, node) )
            vsc.addError("Dynamic config options are only allowed in the `process` scope", node);
        if( inClosure ) {
            vsc.pushScope(ScriptDsl.class);
            if( inProcess )
                vsc.pushScope(ProcessDsl.class);
        }

        super.visitConfigAssign(node);

        if( inClosure ) {
            if( inProcess )
                vsc.popScope();
            vsc.popScope();
        }

        for( int i = 0; i < node.names.size() - 1; i++ )
            configScopes.pop();
    }

    private static boolean isWorkflowHandler(List<String> scopes, ConfigAssignNode node) {
        var option = node.names.get(node.names.size() - 1);
        return scopes.size() == 1
            && "workflow".equals(scopes.get(0))
            && List.of("onComplete", "onError").contains(option);
    }

    @Override
    public void visitConfigBlock(ConfigBlockNode node) {
        var newScope = node.kind == null;
        if( newScope )
            configScopes.add(node.name);
        super.visitConfigBlock(node);
        if( newScope )
            configScopes.remove(configScopes.size() - 1);
    }

    @Override
    public void visitConfigInclude(ConfigIncludeNode node) {
        checkConfigInclude(node);
        visit(node.source);
    }

    private void checkConfigInclude(ConfigIncludeNode node) {
        if( configScopes.isEmpty() )
            return;
        if( configScopes.size() == 2 && "profiles".equals(configScopes.get(0)) )
            return;
        vsc.addError("Config includes are only allowed at the top-level or in a profile", node);
    }

    // statements

    @Override
    public void visitBlockStatement(BlockStatement node) {
        var newScope = node.getVariableScope() != null;
        if( newScope ) vsc.pushScope();
        node.setVariableScope(currentScope());
        super.visitBlockStatement(node);
        if( newScope ) vsc.popScope();
    }

    @Override
    public void visitCatchStatement(CatchStatement node) {
        vsc.pushScope();
        vsc.declare(node.getVariable(), node);
        super.visitCatchStatement(node);
        vsc.popScope();
    }

    @Override
    public void visitExpressionStatement(ExpressionStatement node) {
        var exp = node.getExpression();
        if( exp instanceof BinaryExpression be && Types.isAssignment(be.getOperation().getType()) ) {
            var source = be.getRightExpression();
            var target = be.getLeftExpression();
            visit(source);
            if( !checkImplicitDeclaration(target) ) {
                visit(target);
            }
            return;
        }
        super.visitExpressionStatement(node);
    }

    private boolean checkImplicitDeclaration(Expression node) {
        if( node instanceof TupleExpression te ) {
            var result = false;
            for( var el : te.getExpressions() )
                result |= declareAssignedVariable((VariableExpression) el);
            return result;
        }
        else if( node instanceof VariableExpression ve ) {
            return declareAssignedVariable(ve);
        }
        return false;
    }

    private boolean declareAssignedVariable(VariableExpression ve) {
        var variable = vsc.findVariableDeclaration(ve.getName(), ve);
        if( variable != null ) {
            ve.setAccessedVariable(variable);
            return false;
        }
        else {
            vsc.addError("`" + ve.getName() + "` was assigned but not declared", ve);
            return true;
        }
    }

    // expressions

    private static final List<String> KEYWORDS = List.of(
        "case",
        "for",
        "switch",
        "while"
    );

    @Override
    public void visitMethodCallExpression(MethodCallExpression node) {
        checkMethodCall(node);
        super.visitMethodCallExpression(node);
    }

    private void checkMethodCall(MethodCallExpression node) {
        if( !node.isImplicitThis() )
            return;
        var name = node.getMethodAsString();
        var defNode = vsc.findDslFunction(name, node);
        if( defNode != null )
            node.putNodeMetaData(ASTNodeMarker.METHOD_TARGET, defNode);
        else if( !KEYWORDS.contains(name) )
            vsc.addError("`" + name + "` is not defined", node.getMethod());
    }

    @Override
    public void visitDeclarationExpression(DeclarationExpression node) {
        visit(node.getRightExpression());

        if( node.isMultipleAssignmentDeclaration() ) {
            for( var el : node.getTupleExpression() )
                vsc.declare((VariableExpression) el);
        }
        else {
            vsc.declare(node.getVariableExpression());
        }
    }

    @Override
    public void visitClosureExpression(ClosureExpression node) {
        vsc.pushScope();
        node.setVariableScope(currentScope());
        if( node.getParameters() != null ) {
            for( var parameter : node.getParameters() ) {
                vsc.declare(parameter, parameter);
                if( parameter.hasInitialExpression() )
                    visit(parameter.getInitialExpression());
            }
        }
        super.visitClosureExpression(node);
        for( var it = currentScope().getReferencedLocalVariablesIterator(); it.hasNext(); ) {
            var variable = it.next();
            variable.setClosureSharedVariable(true);
        }
        vsc.popScope();
    }

    @Override
    public void visitVariableExpression(VariableExpression node) {
        var name = node.getName();
        Variable variable = vsc.findVariableDeclaration(name, node);
        if( variable == null ) {
            if( "it".equals(name) ) {
                vsc.addParanoidWarning("Implicit closure parameter `it` will not be supported in a future version", node);
            }
            else {
                variable = new DynamicVariable(name, false);
            }
        }
        if( variable != null ) {
            node.setAccessedVariable(variable);
        }
    }

    // helpers

    private VariableScope currentScope() {
        return vsc.getCurrentScope();
    }

    private List<String> currentConfigScopes() {
        var names = new ArrayList<>(configScopes);
        if( !names.isEmpty() && "profiles".equals(names.get(0)) ) {
            if( !names.isEmpty() ) names.remove(0);
            if( !names.isEmpty() ) names.remove(0);
        }
        return names;
    }

}
