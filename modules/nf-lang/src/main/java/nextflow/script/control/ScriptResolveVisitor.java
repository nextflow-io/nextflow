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
import java.util.Collections;
import java.util.List;

import nextflow.script.ast.AssignmentExpression;
import nextflow.script.ast.FunctionNode;
import nextflow.script.ast.IncludeNode;
import nextflow.script.ast.OutputNode;
import nextflow.script.ast.ParamNodeV1;
import nextflow.script.ast.ProcessNodeV1;
import nextflow.script.ast.ProcessNodeV2;
import nextflow.script.ast.RecordNode;
import nextflow.script.ast.ScriptNode;
import nextflow.script.ast.ScriptVisitorSupport;
import nextflow.script.ast.TupleParameter;
import nextflow.script.ast.WorkflowNode;
import org.codehaus.groovy.ast.ClassNode;
import org.codehaus.groovy.ast.FieldNode;
import org.codehaus.groovy.ast.DynamicVariable;
import org.codehaus.groovy.ast.Parameter;
import org.codehaus.groovy.ast.expr.VariableExpression;
import org.codehaus.groovy.ast.stmt.ExpressionStatement;
import org.codehaus.groovy.ast.stmt.Statement;
import org.codehaus.groovy.control.CompilationUnit;
import org.codehaus.groovy.control.SourceUnit;

import static nextflow.script.ast.ASTUtils.*;

/**
 * Resolve variable names, function names, and type names in
 * a script.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public class ScriptResolveVisitor extends ScriptVisitorSupport {

    private SourceUnit sourceUnit;

    private List<ClassNode> imports;

    private ResolveVisitor resolver;

    public ScriptResolveVisitor(SourceUnit sourceUnit, CompilationUnit compilationUnit, List<ClassNode> defaultImports, List<ClassNode> libImports) {
        this.sourceUnit = sourceUnit;
        this.imports = new ArrayList<>(defaultImports);
        this.resolver = new ResolveVisitor(sourceUnit, compilationUnit, imports, libImports);
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

            // append included types to default imports
            for( var includeNode : sn.getIncludes() ) {
                for( var entry : includeNode.entries ) {
                    if( entry.getTarget() instanceof ClassNode cn )
                        imports.add(cn);
                }
            }
    
            // resolve type names
            if( sn.getParams() != null )
                visitParams(sn.getParams());
            for( var paramNode : sn.getParamsV1() )
                visitParamV1(paramNode);
            for( var workflowNode : sn.getWorkflows() )
                visitWorkflow(workflowNode);
            for( var processNode : sn.getProcesses() )
                visitProcess(processNode);
            for( var functionNode : sn.getFunctions() )
                visitFunction(functionNode);
            for( var type : sn.getClasses() ) {
                if( type instanceof RecordNode rn )
                    visitRecord(rn);
                else if( type.isEnum() )
                    visitEnum(type);
            }
            if( sn.getOutputs() != null )
                visitOutputs(sn.getOutputs());
    
            // report errors for any unresolved variable references
            new DynamicVariablesVisitor().visit(sn);
        }
    }

    @Override
    public void visitParam(Parameter node) {
        node.setInitialExpression(resolver.transform(node.getInitialExpression()));
        resolver.resolveOrFail(node.getType(), node);
    }

    @Override
    public void visitParamV1(ParamNodeV1 node) {
        node.value = resolver.transform(node.value);
    }

    @Override
    public void visitWorkflow(WorkflowNode node) {
        for( var take : node.getParameters() )
            resolver.resolveOrFail(take.getType(), take);
        resolver.visit(node.main);
        resolveTypedOutputs(node.emits);
        resolver.visit(node.emits);
        resolver.visit(node.publishers);
        resolver.visit(node.onComplete);
        resolver.visit(node.onError);
    }

    private void resolveTypedOutputs(Statement block) {
        for( var stmt : asBlockStatements(block) ) {
            var stmtX = (ExpressionStatement)stmt;
            var output = stmtX.getExpression();
            var target =
                output instanceof AssignmentExpression ae ? ae.getLeftExpression() :
                output instanceof VariableExpression ve ? ve :
                null;

            if( target instanceof VariableExpression ve )
                resolver.resolveOrFail(ve);
        }
    }

    @Override
    public void visitProcessV2(ProcessNodeV2 node) {
        for( var input : node.inputs ) {
            resolver.resolveOrFail(input.getType(), input);
            if( input instanceof TupleParameter tp && tp.isRecord() ) {
                for( var component : tp.components )
                    resolver.resolveOrFail(component.getType(), component);
            }
        }
        resolver.visit(node.directives);
        resolver.visit(node.stagers);
        resolveTypedOutputs(node.outputs);
        resolver.visit(node.outputs);
        resolver.visit(node.topics);
        resolver.visit(node.when);
        resolver.visit(node.exec);
        resolver.visit(node.stub);
    }

    @Override
    public void visitProcessV1(ProcessNodeV1 node) {
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
    public void visitField(FieldNode node) {
        resolver.resolveOrFail(node.getType(), node);
    }

    @Override
    public void visitOutput(OutputNode node) {
        resolver.resolveOrFail(node.getType(), node.getType());
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
