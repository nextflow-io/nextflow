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
package nextflow.script.ast;

import org.codehaus.groovy.ast.ClassCodeVisitorSupport;
import org.codehaus.groovy.ast.ClassNode;
import org.codehaus.groovy.ast.Parameter;
import org.codehaus.groovy.ast.expr.ElvisOperatorExpression;
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
        if( script.getParams() != null )
            visitParams(script.getParams());
        for( var paramNode : script.getParamsV1() )
            visitParamV1(paramNode);
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
        if( script.getOutputs() != null )
            visitOutputs(script.getOutputs());
    }

    @Override
    public void visitFeatureFlag(FeatureFlagNode node) {
    }

    @Override
    public void visitInclude(IncludeNode node) {
        visit(node.source);
    }

    @Override
    public void visitParams(ParamBlockNode node) {
        for( var param : node.declarations )
            visitParam(param);
    }

    @Override
    public void visitParam(Parameter node) {
    }

    @Override
    public void visitParamV1(ParamNodeV1 node) {
        visit(node.value);
    }

    @Override
    public void visitWorkflow(WorkflowNode node) {
        visit(node.main);
        visit(node.emits);
        visit(node.publishers);
        visit(node.onComplete);
        visit(node.onError);
    }

    @Override
    public void visitProcess(ProcessNode node) {
        if( node instanceof ProcessNodeV2 pn )
            visitProcessV2(pn);
        if( node instanceof ProcessNodeV1 pn )
            visitProcessV1(pn);
    }

    @Override
    public void visitProcessV2(ProcessNodeV2 node) {
        visit(node.directives);
        visit(node.stagers);
        visit(node.outputs);
        visit(node.topics);
        visit(node.when);
        visit(node.exec);
        visit(node.stub);
    }

    @Override
    public void visitProcessV1(ProcessNodeV1 node) {
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
    public void visitOutputs(OutputBlockNode node) {
        for( var output : node.declarations )
            visitOutput(output);
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

    @Override
    public void visitShortTernaryExpression(ElvisOperatorExpression node) {
        node.getTrueExpression().visit(this);
        node.getFalseExpression().visit(this);
    }

}
