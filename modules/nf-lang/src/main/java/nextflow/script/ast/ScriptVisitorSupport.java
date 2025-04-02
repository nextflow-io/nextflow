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
package nextflow.script.ast;

import org.codehaus.groovy.ast.ClassCodeVisitorSupport;
import org.codehaus.groovy.ast.ClassNode;
import org.codehaus.groovy.ast.expr.MethodCallExpression;

public abstract class ScriptVisitorSupport extends ClassCodeVisitorSupport implements ScriptVisitor {

    //--------------------------------------------------------------------------
    // script declarations

    @Override
    public void visit(ScriptNode script) {
        for( var featureFlag : script.getFeatureFlags() )
            visitFeatureFlag(featureFlag);
        for( var includeNode : script.getIncludes() )
            visitInclude(includeNode);
        for( var paramNode : script.getParams() )
            visitParam(paramNode);
        for( var workflowNode : script.getWorkflows() )
            visitWorkflow(workflowNode);
        for( var processNode : script.getProcesses() )
            visitProcess(processNode);
        for( var functionNode : script.getFunctions() )
            visitFunction(functionNode);
        for( var classNode : script.getClasses() ) {
            if( classNode.isEnum() )
                visitEnum(classNode);
        }
        if( script.getOutput() != null )
            visitOutput(script.getOutput());
    }

    @Override
    public void visitFeatureFlag(FeatureFlagNode node) {
    }

    @Override
    public void visitInclude(IncludeNode node) {
        visit(node.source);
    }

    @Override
    public void visitParam(ParamNode node) {
        visit(node.value);
    }

    @Override
    public void visitWorkflow(WorkflowNode node) {
        visit(node.takes);
        visit(node.main);
        visit(node.emits);
        visit(node.publishers);
    }

    @Override
    public void visitProcess(ProcessNode node) {
        visit(node.directives);
        visit(node.inputs);
        visit(node.outputs);
        visit(node.when);
        visit(node.exec);
        visit(node.stub);
    }

    @Override
    public void visitFunction(FunctionNode node) {
        visit(node.getCode());
    }

    @Override
    public void visitEnum(ClassNode node) {
        for( var fn : node.getFields() )
            visitField(fn);
    }

    @Override
    public void visitOutput(OutputNode node) {
        visit(node.body);
    }

    //--------------------------------------------------------------------------
    // expressions

    @Override
    public void visitMethodCallExpression(MethodCallExpression node) {
        if( !node.isImplicitThis() )
            node.getObjectExpression().visit(this);
        node.getMethod().visit(this);
        node.getArguments().visit(this);
    }

}
