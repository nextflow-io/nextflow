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
package nextflow.script.control;

import java.lang.reflect.Modifier;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import nextflow.script.ast.AgentNode;
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
import nextflow.script.types.Record;
import nextflow.script.types.Tuple;
import org.codehaus.groovy.ast.ClassHelper;
import org.codehaus.groovy.ast.ClassNode;
import org.codehaus.groovy.ast.FieldNode;
import org.codehaus.groovy.ast.DynamicVariable;
import org.codehaus.groovy.ast.GenericsType;
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

    private static final ClassNode RECORD_TYPE = ClassHelper.makeCached(Record.class);
    private static final ClassNode TUPLE_TYPE = ClassHelper.makeCached(Tuple.class);

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
            for( var agentNode : sn.getAgents() )
                visitAgent(agentNode);
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

            // check agent outputs once record field types are resolved
            for( var agentNode : sn.getAgents() )
                checkAgentOutputs(agentNode.outputs);

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

    @Override
    public void visitAgent(AgentNode node) {
        resolveInputs(node.inputs);
        resolver.visit(node.directives);
        resolveTypedOutputs(node.outputs);
        resolver.visit(node.outputs);
        resolver.visit(node.prompt);
    }

    /**
     * A typed agent output (`name: Type`) is answered by the model under
     * a JSON schema, so it must declare a supported type.
     */
    private void checkAgentOutputs(Statement block) {
        for( var stmt : asBlockStatements(block) ) {
            if( !(((ExpressionStatement) stmt).getExpression() instanceof VariableExpression ve) )
                continue;
            // a bare name without a type refers to an existing variable, such as an input
            var av = ve.getAccessedVariable();
            var type = ve.getOriginType();
            if( (av != null && av != ve) || ClassHelper.isDynamicTyped(type) ) {
                resolver.addError("Agent output `" + ve.getName() + "` should declare a type -- typed outputs are answered by the model", ve);
                continue;
            }
            var unsupported = unsupportedAgentOutputType(type, new HashSet<>());
            if( unsupported != null )
                resolver.addError("Agent output `" + ve.getName() + "` has unsupported " + unsupported + " -- supported types are Boolean, Float, Integer, List<E>, Path, String, or a record type", ve);
        }
    }

    // the declared type name, since the display name of e.g. Double is Float
    private static String typeName(ClassNode type) {
        var gts = type.getGenericsTypes();
        if( gts == null )
            return type.getNameWithoutPackage();
        var args = Arrays.stream(gts).map(gt -> typeName(gt.getType())).toList();
        return type.getNameWithoutPackage() + "<" + String.join(", ", args) + ">";
    }

    private static final List<ClassNode> AGENT_OUTPUT_TYPES = List.of(
        ClassHelper.Boolean_TYPE,
        ClassHelper.Float_TYPE,
        ClassHelper.Integer_TYPE,
        ClassHelper.makeCached(java.nio.file.Path.class),
        ClassHelper.STRING_TYPE
    );

    /**
     * Get a description of the unsupported part of an agent output type,
     * or null if the type is supported.
     *
     * @param type
     * @param visited
     */
    private static String unsupportedAgentOutputType(ClassNode type, Set<ClassNode> visited) {
        if( ClassHelper.LIST_TYPE.equals(type) ) {
            var gts = type.getGenericsTypes();
            if( gts == null || gts.length != 1 )
                return "type " + typeName(type);
            var result = unsupportedAgentOutputType(gts[0].getType(), visited);
            return result == null || result.startsWith("field ") ? result : "type " + typeName(type);
        }
        if( type.redirect() instanceof RecordNode rn ) {
            if( !visited.add(rn) )
                return null;
            for( var fn : rn.getFields() ) {
                var result = unsupportedAgentOutputType(fn.getType(), visited);
                if( result != null )
                    return result.startsWith("field ") ? result : "field `" + fn.getName() + "` with " + result;
            }
            return null;
        }
        return AGENT_OUTPUT_TYPES.contains(type) ? null : "type " + typeName(type);
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
        resolveInputs(node.inputs);
        resolver.visit(node.directives);
        resolver.visit(node.stagers);
        resolveTypedOutputs(node.outputs);
        resolver.visit(node.outputs);
        resolver.visit(node.topics);
        resolver.visit(node.when);
        resolver.visit(node.exec);
        resolver.visit(node.stub);
    }

    private void resolveInputs(Parameter[] inputs) {
        for( var input : asFlatParams(inputs) ) {
            resolver.resolveOrFail(input.getType(), input);
        }
        for( var input : inputs ) {
            var type = input.getType();
            if( input instanceof TupleParameter tp && RECORD_TYPE.equals(type) )
                resolveRecordInput(tp);
            if( input instanceof TupleParameter tp && TUPLE_TYPE.equals(type) )
                resolveTupleInput(tp);
        }
    }

    private void resolveRecordInput(TupleParameter tp) {
        var type = tp.getType();
        for( var param : tp.components ) {
            var fn = new FieldNode(param.getName(), Modifier.PUBLIC, param.getType(), type, null);
            fn.setDeclaringClass(type);
            type.addField(fn);
        }
    }

    private void resolveTupleInput(TupleParameter tp) {
        var genericsTypes = Arrays.stream(tp.components)
            .map(p -> new GenericsType(p.getType()))
            .toArray(GenericsType[]::new);
        tp.getType().setGenericsTypes(genericsTypes);
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
