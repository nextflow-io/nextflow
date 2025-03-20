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

import java.util.Collections;
import java.util.HashMap;
import java.util.IdentityHashMap;
import java.util.Map;
import java.util.Set;

import nextflow.script.ast.ASTNodeMarker;
import nextflow.script.dsl.Constant;
import nextflow.script.dsl.Operator;
import nextflow.script.ast.ProcessNode;
import nextflow.script.ast.WorkflowNode;
import org.codehaus.groovy.ast.ASTNode;
import org.codehaus.groovy.ast.ClassHelper;
import org.codehaus.groovy.ast.ClassNode;
import org.codehaus.groovy.ast.FieldNode;
import org.codehaus.groovy.ast.MethodNode;
import org.codehaus.groovy.ast.Parameter;
import org.codehaus.groovy.ast.PropertyNode;
import org.codehaus.groovy.ast.Variable;
import org.codehaus.groovy.ast.VariableScope;
import org.codehaus.groovy.ast.expr.VariableExpression;
import org.codehaus.groovy.control.SourceUnit;
import org.codehaus.groovy.control.messages.SyntaxErrorMessage;
import org.codehaus.groovy.control.messages.WarningMessage;
import org.codehaus.groovy.syntax.SyntaxException;
import org.codehaus.groovy.syntax.Token;
import org.codehaus.groovy.syntax.Types;

import static nextflow.script.ast.ASTUtils.*;

/**
 * Resolve variable and function names.
 *
 * @see org.codehaus.groovy.classgen.VariableScopeVisitor
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public class VariableScopeChecker {

    private SourceUnit sourceUnit;

    private Map<String,MethodNode> includes = new HashMap<>();

    private VariableScope currentScope;

    private Set<Variable> declaredVariables = Collections.newSetFromMap(new IdentityHashMap<>());

    public VariableScopeChecker(SourceUnit sourceUnit, ClassNode classScope) {
        this.sourceUnit = sourceUnit;
        this.currentScope = new VariableScope();
        this.currentScope.setClassScope(classScope);
    }

    public void setCurrentScope(VariableScope scope) {
        currentScope = scope;
    }

    public VariableScope getCurrentScope() {
        return currentScope;
    }

    public void include(String name, MethodNode variable) {
        includes.put(name, variable);
    }

    public MethodNode getInclude(String name) {
        return includes.get(name);
    }

    public void checkUnusedVariables() {
        for( var variable : declaredVariables ) {
            if( variable instanceof ASTNode node && !variable.getName().startsWith("_") ) {
                var message = variable instanceof Parameter
                    ? "Parameter was not used -- prefix with `_` to suppress warning"
                    : "Variable was declared but not used";
                sourceUnit.addWarning(message, node);
            }
        }
    }

    public void pushScope(Class classScope) {
        currentScope = new VariableScope(currentScope);
        if( classScope != null )
            currentScope.setClassScope(ClassHelper.makeCached(classScope));
    }

    public void pushScope() {
        pushScope(null);
    }

    public void popScope() {
        currentScope = currentScope.getParent();
    }

    public void declare(VariableExpression variable) {
        declare(variable, variable);
        variable.setAccessedVariable(variable);
    }

    public void declare(Variable variable, ASTNode context) {
        var name = variable.getName();
        for( var scope = currentScope; scope != null; scope = scope.getParent() ) {
            var other = scope.getDeclaredVariable(name);
            if( other != null ) {
                addError("`" + name + "` is already declared", context, "First declared here", (ASTNode) other);
                break;
            }
        }
        currentScope.putDeclaredVariable(variable);
        declaredVariables.add(variable);
    }

    /**
     * Find the declaration of a given variable.
     *
     * @param name
     * @param node
     */
    public Variable findVariableDeclaration(String name, ASTNode node) {
        Variable variable = null;
        VariableScope scope = currentScope;
        boolean isClassVariable = false;
        while( scope != null ) {
            variable = scope.getDeclaredVariable(name);
            if( variable != null )
                break;
            variable = scope.getReferencedLocalVariable(name);
            if( variable != null )
                break;
            variable = scope.getReferencedClassVariable(name);
            if( variable != null ) {
                isClassVariable = true;
                break;
            }
            variable = findDslVariable(scope.getClassScope(), name, node);
            if( variable != null ) {
                isClassVariable = true;
                break;
            }
            scope = scope.getParent();
        }
        if( variable == null )
            return null;
        VariableScope end = scope;
        scope = currentScope;
        while( true ) {
            if( isClassVariable )
                scope.putReferencedClassVariable(variable);
            else
                scope.putReferencedLocalVariable(variable);
            if( scope == end )
                break;
            scope = scope.getParent();
        }
        declaredVariables.remove(variable);
        return variable;
    }

    /**
     * Find the definition of a built-in variable.
     *
     * @param cn
     * @param name
     * @param node
     */
    private Variable findDslVariable(ClassNode cn, String name, ASTNode node) {
        while( cn != null ) {
            for( var mn : cn.getMethods() ) {
                // processes, workflows, and operators can be accessed as variables, e.g. with pipes
                if( isDataflowMethod(mn) && name.equals(mn.getName()) ) {
                    return wrapMethodAsVariable(mn, name);
                }
                // built-in variables are methods annotated as @Constant
                var an = findAnnotation(mn, Constant.class);
                if( !an.isPresent() )
                    continue;
                if( !name.equals(an.get().getMember("value").getText()) )
                    continue;
                if( findAnnotation(mn, Deprecated.class).isPresent() )
                    addParanoidWarning("`" + name + "` is deprecated and will be removed in a future version", node);
                return wrapMethodAsVariable(mn, name);
            }

            cn = cn.getInterfaces().length > 0
                ? cn.getInterfaces()[0]
                : null;
        }

        if( includes.containsKey(name) )
            return wrapMethodAsVariable(includes.get(name), name);

        return null;
    }

    public static boolean isDataflowMethod(MethodNode mn) {
        return mn instanceof ProcessNode || mn instanceof WorkflowNode || isOperator(mn);
    }

    public static boolean isOperator(MethodNode mn) {
        return findAnnotation(mn, Operator.class).isPresent();
    }

    private static PropertyNode wrapMethodAsVariable(MethodNode mn, String name) {
        var cn = mn.getDeclaringClass();
        var fn = new FieldNode(name, mn.getModifiers() & 0xF, mn.getReturnType(), cn, null);
        fn.setHasNoRealSourcePosition(true);
        fn.setDeclaringClass(cn);
        fn.setSynthetic(true);
        var pn = new PropertyNode(fn, fn.getModifiers(), null, null);
        pn.putNodeMetaData(ASTNodeMarker.METHOD_VARIABLE_TARGET, mn);
        pn.setDeclaringClass(cn);
        return pn;
    }

    /**
     * Find the definition of a built-in function.
     *
     * @param name
     * @param node
     */
    public MethodNode findDslFunction(String name, ASTNode node) {
        VariableScope scope = currentScope;
        while( scope != null ) {
            ClassNode cn = scope.getClassScope();
            while( cn != null ) {
                for( var mn : cn.getMethods() ) {
                    // built-in functions are methods not annotated as @Constant
                    if( findAnnotation(mn, Constant.class).isPresent() )
                        continue;
                    if( !name.equals(mn.getName()) )
                        continue;
                    if( findAnnotation(mn, Deprecated.class).isPresent() )
                        addParanoidWarning("`" + name + "` is deprecated and will be removed in a future version", node);
                    return mn;
                }
    
                cn = cn.getInterfaces().length > 0
                    ? cn.getInterfaces()[0]
                    : null;
            }
            scope = scope.getParent();
        }

        return includes.get(name);
    }

    public void addParanoidWarning(String message, String tokenText, ASTNode node, String otherMessage, ASTNode otherNode) {
        var token = new Token(0, tokenText, node.getLineNumber(), node.getColumnNumber()); // ASTNode to CSTNode
        var warning = new ParanoidWarning(WarningMessage.POSSIBLE_ERRORS, message, token, sourceUnit);
        if( otherNode != null )
            warning.setRelatedInformation(otherMessage, otherNode);
        sourceUnit.getErrorCollector().addWarning(warning);
    }

    public void addParanoidWarning(String message, ASTNode node, String otherMessage, ASTNode otherNode) {
        addParanoidWarning(message, "", node, otherMessage, otherNode);
    }

    public void addParanoidWarning(String message, String tokenText, ASTNode node) {
        addParanoidWarning(message, tokenText, node, null, null);
    }

    public void addParanoidWarning(String message, ASTNode node) {
        addParanoidWarning(message, "", node, null, null);
    }

    public void addError(String message, ASTNode node) {
        addError(new VariableScopeError(message, node));
    }

    public void addError(String message, ASTNode node, String otherMessage, ASTNode otherNode) {
        var cause = new VariableScopeError(message, node);
        if( otherNode != null )
            cause.setRelatedInformation(otherMessage, otherNode);
        addError(cause);
    }

    public void addError(SyntaxException cause) {
        var errorMessage = new SyntaxErrorMessage(cause, sourceUnit);
        sourceUnit.getErrorCollector().addErrorAndContinue(errorMessage);
    }

    private class VariableScopeError extends SyntaxException implements PhaseAware, RelatedInformationAware {

        private String otherMessage;

        private ASTNode otherNode;

        public VariableScopeError(String message, ASTNode node) {
            super(message, node);
        }

        public void setRelatedInformation(String otherMessage, ASTNode otherNode) {
            this.otherMessage = otherMessage;
            this.otherNode = otherNode;
        }

        @Override
        public int getPhase() {
            return Phases.NAME_RESOLUTION;
        }

        @Override
        public String getOtherMessage() {
            return otherMessage;
        }

        @Override
        public ASTNode getOtherNode() {
            return otherNode;
        }
    }

}
