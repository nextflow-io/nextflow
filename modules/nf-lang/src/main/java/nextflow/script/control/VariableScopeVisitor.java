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

import java.util.HashMap;
import java.util.List;

import nextflow.script.ast.ASTNodeMarker;
import nextflow.script.ast.AssignmentExpression;
import nextflow.script.ast.FeatureFlagNode;
import nextflow.script.ast.FunctionNode;
import nextflow.script.ast.IncludeNode;
import nextflow.script.ast.OutputNode;
import nextflow.script.ast.ProcessNode;
import nextflow.script.ast.ScriptNode;
import nextflow.script.ast.ScriptVisitorSupport;
import nextflow.script.ast.WorkflowNode;
import nextflow.script.dsl.Constant;
import nextflow.script.dsl.EntryWorkflowDsl;
import nextflow.script.dsl.FeatureFlag;
import nextflow.script.dsl.FeatureFlagDsl;
import nextflow.script.dsl.OutputDsl;
import nextflow.script.dsl.ProcessDsl;
import nextflow.script.dsl.ScriptDsl;
import nextflow.script.dsl.WorkflowDsl;
import org.codehaus.groovy.ast.ASTNode;
import org.codehaus.groovy.ast.ClassHelper;
import org.codehaus.groovy.ast.ClassNode;
import org.codehaus.groovy.ast.DynamicVariable;
import org.codehaus.groovy.ast.MethodNode;
import org.codehaus.groovy.ast.Variable;
import org.codehaus.groovy.ast.VariableScope;
import org.codehaus.groovy.ast.expr.BinaryExpression;
import org.codehaus.groovy.ast.expr.ClosureExpression;
import org.codehaus.groovy.ast.expr.ConstantExpression;
import org.codehaus.groovy.ast.expr.DeclarationExpression;
import org.codehaus.groovy.ast.expr.EmptyExpression;
import org.codehaus.groovy.ast.expr.Expression;
import org.codehaus.groovy.ast.expr.MapEntryExpression;
import org.codehaus.groovy.ast.expr.MethodCallExpression;
import org.codehaus.groovy.ast.expr.PropertyExpression;
import org.codehaus.groovy.ast.expr.TupleExpression;
import org.codehaus.groovy.ast.expr.VariableExpression;
import org.codehaus.groovy.ast.stmt.BlockStatement;
import org.codehaus.groovy.ast.stmt.CatchStatement;
import org.codehaus.groovy.ast.stmt.ExpressionStatement;
import org.codehaus.groovy.ast.stmt.Statement;
import org.codehaus.groovy.control.SourceUnit;
import org.codehaus.groovy.control.messages.SyntaxErrorMessage;
import org.codehaus.groovy.syntax.SyntaxException;
import org.codehaus.groovy.syntax.Types;

import static nextflow.script.ast.ASTUtils.*;

/**
 * Initialize the variable scopes for an AST.
 *
 * @see org.codehaus.groovy.classgen.VariableScopeVisitor
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class VariableScopeVisitor extends ScriptVisitorSupport {

    private SourceUnit sourceUnit;

    private VariableScopeChecker vsc;

    private MethodNode currentDefinition;

    public VariableScopeVisitor(SourceUnit sourceUnit) {
        this.sourceUnit = sourceUnit;
        this.vsc = new VariableScopeChecker(sourceUnit, new ClassNode(ScriptDsl.class));
    }

    @Override
    protected SourceUnit getSourceUnit() {
        return sourceUnit;
    }

    public void declare() {
        var moduleNode = sourceUnit.getAST();
        if( moduleNode instanceof ScriptNode sn ) {
            for( var includeNode : sn.getIncludes() )
                declareInclude(includeNode);
            for( var workflowNode : sn.getWorkflows() ) {
                if( !workflowNode.isEntry() )
                    declareMethod(workflowNode);
            }
            for( var processNode : sn.getProcesses() )
                declareMethod(processNode);
            for( var functionNode : sn.getFunctions() )
                declareMethod(functionNode);
        }
    }

    private void declareInclude(IncludeNode node) {
        for( var module : node.modules ) {
            if( module.getTarget() == null )
                continue;
            var name = module.getNameOrAlias();
            var otherInclude = vsc.getInclude(name);
            if( otherInclude != null )
                vsc.addError("`" + name + "` is already included", node, "First included here", otherInclude);
            vsc.include(name, module.getTarget());
        }
    }

    private void declareMethod(MethodNode mn) {
        var cn = currentScope().getClassScope();
        var name = mn.getName();
        var otherInclude = vsc.getInclude(name);
        if( otherInclude != null ) {
            vsc.addError("`" + name + "` is already included", mn, "First included here", otherInclude);
        }
        var otherMethods = cn.getDeclaredMethods(name);
        if( otherMethods.size() > 0 ) {
            var other = otherMethods.get(0);
            var first = mn.getLineNumber() < other.getLineNumber() ? mn : other;
            var second = mn.getLineNumber() < other.getLineNumber() ? other : mn;
            vsc.addError("`" + name + "` is already declared", second, "First declared here", first);
            return;
        }
        cn.addMethod(mn);
    }

    public void visit() {
        var moduleNode = sourceUnit.getAST();
        if( moduleNode instanceof ScriptNode sn ) {
            super.visit(sn);
            vsc.checkUnusedVariables();
        }
    }

    @Override
    public void visitFeatureFlag(FeatureFlagNode node) {
        var cn = ClassHelper.makeCached(FeatureFlagDsl.class);
        var result = cn.getFields().stream()
            .filter(fn ->
                findAnnotation(fn, FeatureFlag.class)
                    .map(an -> an.getMember("value").getText())
                    .map(name -> name.equals(node.name))
                    .orElse(false)
            )
            .findFirst();

        if( result.isPresent() ) {
            var ffn = result.get();
            if( findAnnotation(ffn, Deprecated.class).isPresent() )
                vsc.addParanoidWarning("`" + node.name + "` is deprecated and will be removed in a future version", node.name, node);
            node.target = ffn;
        }
        else {
            vsc.addError("Unrecognized feature flag '" + node.name + "'", node);
        }
    }

    private boolean inWorkflowEmit;

    @Override
    public void visitWorkflow(WorkflowNode node) {
        vsc.pushScope(node.isEntry() ? EntryWorkflowDsl.class : WorkflowDsl.class);
        currentDefinition = node;
        node.setVariableScope(currentScope());

        declareWorkflowInputs(node.takes);

        visit(node.main);
        if( node.main instanceof BlockStatement block )
            copyVariableScope(block.getVariableScope());

        visitWorkflowEmits(node.emits);
        visit(node.publishers);

        currentDefinition = null;
        vsc.popScope();
    }

    private void declareWorkflowInputs(Statement takes) {
        for( var stmt : asBlockStatements(takes) ) {
            var ve = asVarX(stmt);
            if( ve == null )
                continue;
            vsc.declare(ve);
        }
    }

    private void copyVariableScope(VariableScope source) {
        for( var it = source.getDeclaredVariablesIterator(); it.hasNext(); ) {
            var variable = it.next();
            currentScope().putDeclaredVariable(variable);
        }
    }

    private void visitWorkflowEmits(Statement emits) {
        var declaredEmits = new HashMap<String,ASTNode>();
        for( var stmt : asBlockStatements(emits) ) {
            var stmtX = (ExpressionStatement)stmt;
            var emit = stmtX.getExpression();
            if( emit instanceof AssignmentExpression assign ) {
                visit(assign.getRightExpression());

                var target = (VariableExpression)assign.getLeftExpression();
                var name = target.getName();
                var other = declaredEmits.get(name);
                if( other != null )
                    vsc.addError("Workflow emit `" + name + "` is already declared", target, "First declared here", other);
                else
                    declaredEmits.put(name, target);
            }
            else {
                visit(emit);
            }
        }
    }

    @Override
    public void visitProcess(ProcessNode node) {
        vsc.pushScope(ProcessDsl.class);
        currentDefinition = node;
        node.setVariableScope(currentScope());

        declareProcessInputs(node.inputs);

        vsc.pushScope(ProcessDsl.InputDsl.class);
        visitDirectives(node.inputs, "process input qualifier", false);
        vsc.popScope();

        if( !(node.when instanceof EmptyExpression) )
            vsc.addParanoidWarning("Process `when` section will not be supported in a future version", node.when);
        visit(node.when);

        visit(node.exec);
        visit(node.stub);

        vsc.pushScope(ProcessDsl.DirectiveDsl.class);
        visitDirectives(node.directives, "process directive", false);
        vsc.popScope();

        vsc.pushScope(ProcessDsl.OutputDsl.class);
        visitDirectives(node.outputs, "process output qualifier", false);
        vsc.popScope();

        currentDefinition = null;
        vsc.popScope();
    }

    private void declareProcessInputs(Statement inputs) {
        for( var stmt : asBlockStatements(inputs) ) {
            var call = asMethodCallX(stmt);
            if( call == null )
                continue;
            if( "tuple".equals(call.getMethodAsString()) ) {
                for( var arg : asMethodCallArguments(call) ) {
                    if( arg instanceof MethodCallExpression mce )
                        declareProcessInput(mce);
                }
            }
            else if( "each".equals(call.getMethodAsString()) ) {
                var args = asMethodCallArguments(call);
                if( args.size() != 1 )
                    continue;
                var firstArg = args.get(0);
                if( firstArg instanceof MethodCallExpression mce )
                    declareProcessInput(mce);
                else if( firstArg instanceof VariableExpression ve )
                    vsc.declare(ve);
            }
            else {
                declareProcessInput(call);
            }
        }
    }

    private static final List<String> DECLARING_INPUT_TYPES = List.of("val", "file", "path");

    private void declareProcessInput(MethodCallExpression call) {
        if( !DECLARING_INPUT_TYPES.contains(call.getMethodAsString()) )
            return;
        var args = asMethodCallArguments(call);
        if( args.isEmpty() )
            return;
        if( args.get(args.size() - 1) instanceof VariableExpression ve )
            vsc.declare(ve);
    }

    private void visitDirectives(Statement node, String typeLabel, boolean checkSyntaxErrors) {
        if( node instanceof BlockStatement block )
            block.setVariableScope(currentScope());
        for( var stmt : asBlockStatements(node) ) {
            var call = checkDirective(stmt, typeLabel, checkSyntaxErrors);
            if( call != null )
                super.visitMethodCallExpression(call);
        }
    }

    private MethodCallExpression checkDirective(Statement node, String typeLabel, boolean checkSyntaxErrors) {
        var call = asMethodCallX(node);
        if( call == null ) {
            if( checkSyntaxErrors )
                addSyntaxError("Invalid " + typeLabel, node);
            return null;
        }
        var name = call.getMethodAsString();
        var mn = vsc.findDslFunction(name, call.getMethod());
        if( mn != null )
            call.putNodeMetaData(ASTNodeMarker.METHOD_TARGET, mn);
        else
            vsc.addError("Invalid " + typeLabel + " `" + name + "`", node);
        return call;
    }

    private static final List<String> EMIT_AND_TOPIC = List.of("emit", "topic");

    @Override
    public void visitMapEntryExpression(MapEntryExpression node) {
        var classScope = currentScope().getClassScope();
        if( classScope != null && classScope.getTypeClass() == ProcessDsl.OutputDsl.class ) {
            var key = node.getKeyExpression();
            if( key instanceof ConstantExpression && EMIT_AND_TOPIC.contains(key.getText()) )
                return;
        }
        super.visitMapEntryExpression(node);
    }

    @Override
    public void visitFunction(FunctionNode node) {
        vsc.pushScope();
        currentDefinition = node;
        node.setVariableScope(currentScope());

        for( var parameter : node.getParameters() ) {
            if( parameter.hasInitialExpression() )
                visit(parameter.getInitialExpression());
            vsc.declare(parameter, node);
        }
        visit(node.getCode());

        currentDefinition = null;
        vsc.popScope();
    }

    @Override
    public void visitOutput(OutputNode node) {
        if( node.body instanceof BlockStatement block )
            visitOutputBody(block);
    }

    private void visitOutputBody(BlockStatement block) {
        block.setVariableScope(currentScope());

        asDirectives(block).forEach((call) -> {
            var code = asDslBlock(call, 1);
            if( code != null )
                visitTargetBody(code);
        });
    }

    private void visitTargetBody(BlockStatement block) {
        vsc.pushScope(OutputDsl.class);
        block.setVariableScope(currentScope());

        asBlockStatements(block).forEach((stmt) -> {
            // validate target directive
            var call = checkDirective(stmt, "output target directive", true);
            if( call == null )
                return;

            // treat as index definition
            var name = call.getMethodAsString();
            if( "index".equals(name) ) {
                var code = asDslBlock(call, 1);
                if( code != null ) {
                    vsc.pushScope(OutputDsl.IndexDsl.class);
                    visitDirectives(code, "output index directive", true);
                    vsc.popScope();
                    return;
                }
            }

            // treat as regular directive
            super.visitMethodCallExpression(call);
        });
        vsc.popScope();
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
        if( exp instanceof AssignmentExpression ae ) {
            var source = ae.getRightExpression();
            var target = ae.getLeftExpression();
            visit(source);
            if( checkImplicitDeclaration(target) ) {
                ae.putNodeMetaData(ASTNodeMarker.IMPLICIT_DECLARATION, Boolean.TRUE);
            }
            else {
                visitMutatedVariable(target);
                visit(target);
            }
            return;
        }
        super.visitExpressionStatement(node);
    }

    /**
     * In processes and workflows, variables can be declared without `def`
     * and are treated as variables scoped to the process or workflow.
     *
     * @param target
     */
    private boolean checkImplicitDeclaration(Expression target) {
        if( target instanceof TupleExpression te ) {
            var result = false;
            for( var el : te.getExpressions() )
                result |= declareAssignedVariable((VariableExpression) el);
            return result;
        }
        else if( target instanceof VariableExpression ve ) {
            return declareAssignedVariable(ve);
        }
        return false;
    }

    private boolean declareAssignedVariable(VariableExpression ve) {
        var variable = vsc.findVariableDeclaration(ve.getName(), ve);
        if( variable != null ) {
            if( isDslVariable(variable) )
                vsc.addError("Built-in variable cannot be re-assigned", ve);
            ve.setAccessedVariable(variable);
            return false;
        }
        else if( currentDefinition instanceof ProcessNode || currentDefinition instanceof WorkflowNode ) {
            if( currentClosure != null )
                addError("Variables in a closure should be declared with `def`", ve);
            var scope = currentScope();
            currentScope(currentDefinition.getVariableScope());
            vsc.declare(ve);
            currentScope(scope);
            return true;
        }
        else {
            vsc.addError("`" + ve.getName() + "` was assigned but not declared", ve);
            return true;
        }
    }

    private boolean isDslVariable(Variable variable) {
        var mn = asMethodVariable(variable);
        return mn != null && findAnnotation(mn, Constant.class).isPresent();
    }

    private void visitMutatedVariable(Expression node) {
        VariableExpression target = null;
        while( true ) {
            // e.g. obj.prop = 123
            if( node instanceof PropertyExpression pe ) {
                node = pe.getObjectExpression();
            }
            // e.g. list[1] = 123 OR map['a'] = 123
            else if( node instanceof BinaryExpression be && be.getOperation().getType() == Types.LEFT_SQUARE_BRACKET ) {
                node = be.getLeftExpression();
            }
            else {
                if( node instanceof VariableExpression ve )
                    target = ve;
                break;
            }
        }
        if( target == null )
            return;
        var variable = vsc.findVariableDeclaration(target.getName(), target);
        if( isDslVariable(variable) ) {
            if( "params".equals(variable.getName()) )
                sourceUnit.addWarning("Params should be declared at the top-level (i.e. outside the workflow)", target);
            // TODO: re-enable after workflow.onComplete bug is fixed
            // else
            //     addError("Built-in variable cannot be mutated", target);
        }
        else if( variable != null ) {
            checkExternalWriteInAsyncClosure(target, variable);
        }
    }

    private void checkExternalWriteInAsyncClosure(VariableExpression target, Variable variable) {
        if( !(currentDefinition instanceof WorkflowNode) )
            return;
        if( currentClosure == null )
            return;
        var scope = currentClosure.getVariableScope();
        var name = variable.getName();
        if( inOperatorCall && scope.isReferencedLocalVariable(name) && scope.getDeclaredVariable(name) == null )
            sourceUnit.addWarning("Mutating an external variable in an operator closure can lead to a race condition", target);
    }

    // expressions

    private static final List<String> KEYWORDS = List.of(
        "case",
        "for",
        "switch",
        "while"
    );

    private boolean inOperatorCall;

    @Override
    public void visitMethodCallExpression(MethodCallExpression node) {
        var target = checkSetAssignment(node);
        if( target != null ) {
            visit(node.getObjectExpression());
            declareAssignedVariable(target);
            return;
        }
        checkMethodCall(node);
        var ioc = inOperatorCall;
        inOperatorCall = isOperatorCall(node);
        super.visitMethodCallExpression(node);
        inOperatorCall = ioc;
    }

    private static boolean isOperatorCall(MethodCallExpression node) {
        return node.getNodeMetaData(ASTNodeMarker.METHOD_TARGET) instanceof MethodNode mn
            && VariableScopeChecker.isOperator(mn);
    }

    /**
     * Treat `set` and `tap` operators as assignments.
     */
    private VariableExpression checkSetAssignment(MethodCallExpression node) {
        if( !(currentDefinition instanceof WorkflowNode) )
            return null;
        var name = node.getMethodAsString();
        if( !"set".equals(name) && !"tap".equals(name) )
            return null;
        var code = asDslBlock(node, 1);
        if( code == null || code.getStatements().size() != 1 )
            return null;
        return asVarX(code.getStatements().get(0));
    }

    private void checkMethodCall(MethodCallExpression node) {
        if( !node.isImplicitThis() )
            return;
        var name = node.getMethodAsString();
        var mn = vsc.findDslFunction(name, node.getMethod());
        if( mn != null ) {
            if( VariableScopeChecker.isDataflowMethod(mn) )
                checkDataflowMethod(node, mn);
            node.putNodeMetaData(ASTNodeMarker.METHOD_TARGET, mn);
        }
        else if( !KEYWORDS.contains(name) ) {
            vsc.addError("`" + name + "` is not defined", node.getMethod());
        }
    }

    private void checkDataflowMethod(MethodCallExpression node, MethodNode mn) {
        if( !(currentDefinition instanceof WorkflowNode) ) {
            var type = dataflowMethodType(mn);
            vsc.addError(type + " can only be called from a workflow", node);
            return;
        }
        if( currentClosure != null ) {
            var type = dataflowMethodType(mn);
            vsc.addError(type + " cannot be called from within a closure", node);
            return;
        }
    }

    private static String dataflowMethodType(MethodNode mn) {
        if( mn instanceof ProcessNode )
            return "Processes";
        if( mn instanceof WorkflowNode )
            return "Workflows";
        return "Operators";
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

    private ClosureExpression currentClosure;

    @Override
    public void visitClosureExpression(ClosureExpression node) {
        var cl = currentClosure;
        currentClosure = node;

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

        currentClosure = cl;
    }

    @Override
    public void visitVariableExpression(VariableExpression node) {
        var name = node.getName();
        Variable variable = vsc.findVariableDeclaration(name, node);
        if( variable == null ) {
            if( "it".equals(name) ) {
                vsc.addParanoidWarning("Implicit closure parameter `it` will not be supported in a future version", node);
            }
            else if( "args".equals(name) ) {
                vsc.addParanoidWarning("The use of `args` outside the entry workflow will not be supported in a future version", node);
            }
            else if( "params".equals(name) ) {
                vsc.addParanoidWarning("The use of `params` outside the entry workflow will not be supported in a future version", node);
            }
            else {
                variable = new DynamicVariable(name, false);
            }
        }
        if( variable != null ) {
            checkGlobalVariableInProcess(variable, node);
            node.setAccessedVariable(variable);
        }
    }

    private static final List<String> WARN_GLOBALS = List.of(
        "baseDir",
        "launchDir",
        "projectDir",
        "workDir"
    );

    private void checkGlobalVariableInProcess(Variable variable, ASTNode context) {
        if( !(currentDefinition instanceof ProcessNode) )
            return;
        var mn = asMethodVariable(variable);
        if( mn != null && mn.getDeclaringClass().getTypeClass() == ScriptDsl.class ) {
            if( WARN_GLOBALS.contains(variable.getName()) )
                sourceUnit.addWarning("The use of `" + variable.getName() + "` in a process is discouraged -- input files should be provided as process inputs", context);
        }
    }

    // helpers

    private VariableScope currentScope() {
        return vsc.getCurrentScope();
    }

    private void currentScope(VariableScope scope) {
        vsc.setCurrentScope(scope);
    }

    public void addSyntaxError(String message, ASTNode node) {
        var cause = new SyntaxException(message, node);
        var errorMessage = new SyntaxErrorMessage(cause, sourceUnit);
        sourceUnit.getErrorCollector().addErrorAndContinue(errorMessage);
    }

}
