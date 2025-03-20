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
import java.util.List;

import nextflow.script.ast.FunctionNode;
import nextflow.script.ast.OutputNode;
import nextflow.script.ast.ParamNode;
import nextflow.script.ast.ProcessNode;
import nextflow.script.ast.ScriptNode;
import nextflow.script.ast.ScriptVisitorSupport;
import nextflow.script.ast.WorkflowNode;
import org.codehaus.groovy.ast.ClassNode;
import org.codehaus.groovy.ast.DynamicVariable;
import org.codehaus.groovy.ast.expr.VariableExpression;
import org.codehaus.groovy.control.CompilationUnit;
import org.codehaus.groovy.control.SourceUnit;

/**
 * Resolve variable names, function names, and type names in
 * a script.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public class ScriptResolveVisitor extends ScriptVisitorSupport {

    private SourceUnit sourceUnit;

    private ResolveVisitor resolver;

    public ScriptResolveVisitor(SourceUnit sourceUnit, CompilationUnit compilationUnit, List<ClassNode> defaultImports, List<ClassNode> libImports) {
        this.sourceUnit = sourceUnit;
        this.resolver = new ResolveVisitor(sourceUnit, compilationUnit, defaultImports, libImports);
    }

    @Override
    protected SourceUnit getSourceUnit() {
        return sourceUnit;
    }

    public void visit() {
        var moduleNode = sourceUnit.getAST();
        if( moduleNode instanceof ScriptNode sn ) {
            // initialize variable scopes
            var variableScopeVisitor = new VariableScopeVisitor(sourceUnit);
            variableScopeVisitor.declare();
            variableScopeVisitor.visit();
    
            // resolve type names
            for( var paramNode : sn.getParams() )
                visitParam(paramNode);
            for( var workflowNode : sn.getWorkflows() )
                visitWorkflow(workflowNode);
            for( var processNode : sn.getProcesses() )
                visitProcess(processNode);
            for( var functionNode : sn.getFunctions() )
                visitFunction(functionNode);
            if( sn.getOutput() != null )
                visitOutput(sn.getOutput());
    
            // report errors for any unresolved variable references
            new DynamicVariablesVisitor().visit(sn);
        }
    }

    @Override
    public void visitParam(ParamNode node) {
        node.value = resolver.transform(node.value);
    }

    @Override
    public void visitWorkflow(WorkflowNode node) {
        resolver.visit(node.main);
        resolver.visit(node.emits);
        resolver.visit(node.publishers);
    }

    @Override
    public void visitProcess(ProcessNode node) {
        resolver.visit(node.directives);
        resolver.visit(node.inputs);
        resolver.visit(node.outputs);
        resolver.visit(node.when);
        resolver.visit(node.exec);
        resolver.visit(node.stub);
    }

    @Override
    public void visitFunction(FunctionNode node) {
        for( var param : node.getParameters() ) {
            param.setInitialExpression(resolver.transform(param.getInitialExpression()));
            resolver.resolveOrFail(param.getType(), param.getType());
        }
        resolver.resolveOrFail(node.getReturnType(), node);
        resolver.visit(node.getCode());
    }

    @Override
    public void visitOutput(OutputNode node) {
        resolver.visit(node.body);
    }

    private class DynamicVariablesVisitor extends ScriptVisitorSupport {

        @Override
        protected SourceUnit getSourceUnit() {
            return sourceUnit;
        }

        @Override
        public void visitVariableExpression(VariableExpression node) {
            var variable = node.getAccessedVariable();
            if( variable instanceof DynamicVariable )
                resolver.addError("`" + node.getName() + "` is not defined", node);
        }
    }

}
