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

import java.net.URI;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import nextflow.module.spi.RemoteModuleResolverProvider;
import nextflow.script.ast.FunctionNode;
import nextflow.script.ast.IncludeEntryNode;
import nextflow.script.ast.ProcessNodeV1;
import nextflow.script.ast.ProcessNodeV2;
import nextflow.script.ast.IncludeNode;
import nextflow.script.ast.ScriptNode;
import nextflow.script.ast.RecordNode;
import nextflow.script.ast.ScriptVisitorSupport;
import nextflow.script.ast.WorkflowNode;
import nextflow.script.dsl.Nullable;
import org.codehaus.groovy.ast.ASTNode;
import org.codehaus.groovy.ast.AnnotatedNode;
import org.codehaus.groovy.ast.AnnotationNode;
import org.codehaus.groovy.ast.ClassHelper;
import org.codehaus.groovy.ast.ClassNode;
import org.codehaus.groovy.ast.CodeVisitorSupport;
import org.codehaus.groovy.ast.FieldNode;
import org.codehaus.groovy.ast.MethodNode;
import org.codehaus.groovy.ast.Parameter;
import org.codehaus.groovy.ast.expr.ClosureExpression;
import org.codehaus.groovy.ast.expr.DeclarationExpression;
import org.codehaus.groovy.ast.expr.Expression;
import org.codehaus.groovy.ast.expr.VariableExpression;
import org.codehaus.groovy.ast.stmt.BlockStatement;
import org.codehaus.groovy.ast.stmt.ExpressionStatement;
import org.codehaus.groovy.ast.stmt.Statement;
import org.codehaus.groovy.control.SourceUnit;
import org.codehaus.groovy.control.messages.SyntaxErrorMessage;
import org.codehaus.groovy.syntax.SyntaxException;

/**
 * Resolve includes against included source files.
 *
 * This visitor should be applied only after all source files
 * have been parsed.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public class ResolveIncludeVisitor extends ScriptVisitorSupport {

    private SourceUnit sourceUnit;

    private URI uri;

    private Path projectDir;

    private Compiler compiler;

    private Set<URI> changedUris;

    private List<SyntaxErrorMessage> errors = new ArrayList<>();

    private boolean changed;

    public ResolveIncludeVisitor(SourceUnit sourceUnit, Path projectDir, Compiler compiler, Set<URI> changedUris) {
        this.sourceUnit = sourceUnit;
        this.uri = sourceUnit.getSource().getURI();
        this.compiler = compiler;
        this.changedUris = changedUris;
        this.projectDir = projectDir;
    }

    public ResolveIncludeVisitor(SourceUnit sourceUnit, Path projectDir, Compiler compiler) {
        this(sourceUnit, projectDir, compiler, null);
    }

    @Override
    protected SourceUnit getSourceUnit() {
        return sourceUnit;
    }

    public void visit() {
        var moduleNode = sourceUnit.getAST();
        if( moduleNode instanceof ScriptNode sn )
            super.visit(sn);
    }

    @Override
    public void visitInclude(IncludeNode node) {
        var source = node.source.getText();
        if( source.startsWith("plugin/") ) {
            setPlaceholderTargets(node);
            return;
        }

        URI includeUri;
        try {
            includeUri = getIncludeUri(uri, source);
        }
        catch( Exception e ) {
            addError(e.getMessage(), node);
            return;
        }

        if( !isIncludeStale(node, includeUri) )
            return;
        changed = true;
        for( var entry : node.entries )
            entry.setTarget(null);
        var includeUnit = compiler.getSource(includeUri);
        if( includeUnit == null ) {
            addError("Invalid include source: '" + includeUri.getPath() + "'", node);
            return;
        }
        if( includeUnit.getAST() == null ) {
            addError("Module could not be parsed: '" + includeUri.getPath() + "'", node);
            return;
        }
        var scriptNode = (ScriptNode) includeUnit.getAST();
        var definitions = getDefinitions(includeUri);
        var pipelineEntries = 0;
        var hasPipeline = false;
        var hasOutput = false;
        for( var entry : node.entries ) {
            var includedName = entry.name;
            var definitionNode = definitions.stream()
                .filter(defNode -> includedName.equals(definitionName(defNode)))
                .findFirst()
                .orElse(null);
            // a `params` or `output` entry that doesn't match a definition of the
            // module refers to the corresponding block of the pipeline
            var blockNode = definitionNode == null ? pipelineBlockType(scriptNode, entry) : null;
            var includedNode = definitionNode != null ? (AnnotatedNode) definitionNode : blockNode;
            if( includedNode == null ) {
                addError("Included name '" + includedName + "' is not defined in module '" + includeUri.getPath() + "'", node);
                continue;
            }
            // a definition that happens to be named `params`, `workflow` or
            // `output` is included as itself, not as a block of the pipeline
            var isPipelineBlock = blockNode != null || isEntryWorkflow(includedName, definitionNode);
            if( isPipelineBlock ) {
                if( entry.alias == null ) {
                    addError("An included pipeline must be aliased, e.g. `" + includedName + " as MY_PIPELINE`", node);
                    continue;
                }
                if( PIPELINE_NAME.equals(includedName) && ++pipelineEntries > 1 ) {
                    addError("A pipeline can be included only once in an include declaration -- use a separate include declaration for each alias", node);
                    continue;
                }
                hasPipeline |= PIPELINE_NAME.equals(includedName);
                hasOutput |= OUTPUT_NAME.equals(includedName);
            }
            entry.setTarget(includedNode);
        }
        // the output block of an included pipeline is resolved against the params
        // of the pipeline, so it can only be included alongside the pipeline itself
        if( hasOutput && !hasPipeline )
            addError("The output block of a pipeline can be included only alongside the pipeline itself -- add `" + PIPELINE_NAME + " as MY_PIPELINE` to this include declaration", node);
        if( hasPipeline )
            checkPipelineParams(node, includeUri);
    }

    private static boolean isEntryWorkflow(String name, AnnotatedNode node) {
        return PIPELINE_NAME.equals(name) && node instanceof WorkflowNode wn && wn.isEntry();
    }

    /**
     * A pipeline receives its params when it is called, so `params` is only
     * available in its entry workflow and output block. A pipeline that refers
     * to `params` anywhere else cannot be included, because that reference
     * would resolve against the calling pipeline instead.
     *
     * The check covers the modules included by the pipeline, because a
     * reference in a module resolves against the calling pipeline just the
     * same, without even a warning.
     */
    private void checkPipelineParams(IncludeNode node, URI includeUri) {
        var offending = findParamsRef(includeUri, new HashSet<>());
        if( offending != null )
            addError("An included pipeline cannot refer to `params` outside of its entry workflow and output block -- declare a workflow or process input instead (in module '" + offending + "')", node);
    }

    /**
     * Find a reference to `params` in a module or the modules that it includes,
     * and return the path of the offending module.
     *
     * @param moduleUri
     * @param visited
     */
    private String findParamsRef(URI moduleUri, Set<URI> visited) {
        if( !visited.add(moduleUri) )
            return null;
        var unit = compiler.getSource(moduleUri);
        if( unit == null || !(unit.getAST() instanceof ScriptNode sn) )
            return null;
        if( hasParamsRef(sn) )
            return moduleUri.getPath();
        for( var include : sn.getIncludes() ) {
            var source = include.source.getText();
            if( source.startsWith("plugin/") )
                continue;
            URI childUri;
            try {
                childUri = getIncludeUri(moduleUri, source);
            }
            catch( Exception e ) {
                continue;
            }
            var result = findParamsRef(childUri, visited);
            if( result != null )
                return result;
        }
        return null;
    }

    private static boolean hasParamsRef(ScriptNode sn) {
        var visitor = new ParamsRefVisitor();
        for( var wn : sn.getWorkflows() ) {
            if( wn.isEntry() )
                continue;
            visitor.visitDeclaration(wn.getParameters(), wn.main, wn.emits, wn.publishers, wn.onComplete, wn.onError);
        }
        for( var pn : sn.getProcesses() ) {
            if( pn instanceof ProcessNodeV1 p1 )
                visitor.visitDeclaration(Parameter.EMPTY_ARRAY, p1.directives, p1.inputs, p1.outputs, asStatement(p1.when), p1.exec, p1.stub);
            else if( pn instanceof ProcessNodeV2 p2 )
                visitor.visitDeclaration(p2.inputs, p2.directives, p2.stagers, p2.outputs, p2.topics, asStatement(p2.when), p2.exec, p2.stub);
        }
        for( var an : sn.getAgents() )
            visitor.visitDeclaration(an.inputs, an.directives, an.outputs, an.prompt);
        for( var fn : sn.getFunctions() )
            visitor.visitDeclaration(fn.getParameters(), fn.getCode());
        return visitor.found;
    }

    private static Statement asStatement(Expression node) {
        return node != null ? new ExpressionStatement(node) : null;
    }

    /**
     * Finds a reference to `params` that is not shadowed by a local binding
     * of the same name.
     */
    private static class ParamsRefVisitor extends CodeVisitorSupport {

        private static final String PARAMS = "params";

        public boolean found;

        private boolean shadowed;

        /**
         * Visit a declaration -- a workflow, process, agent or function -- and
         * the default values of its inputs.
         *
         * An input named `params` shadows the params of the pipeline, so a
         * reference within the declaration body is not an offending reference.
         *
         * @param parameters
         * @param body
         */
        public void visitDeclaration(Parameter[] parameters, Statement... body) {
            shadowed = false;
            for( var parameter : parameters ) {
                if( parameter.hasInitialExpression() )
                    parameter.getInitialExpression().visit(this);
            }
            if( isShadowing(parameters) )
                return;
            for( var node : body ) {
                if( node != null )
                    node.visit(this);
            }
        }

        @Override
        public void visitBlockStatement(BlockStatement node) {
            final var saved = shadowed;
            super.visitBlockStatement(node);
            shadowed = saved;
        }

        @Override
        public void visitClosureExpression(ClosureExpression node) {
            if( isShadowing(node.getParameters()) )
                return;
            final var saved = shadowed;
            super.visitClosureExpression(node);
            shadowed = saved;
        }

        @Override
        public void visitDeclarationExpression(DeclarationExpression node) {
            // a local variable named `params` shadows the pipeline params for
            // the remainder of the block that declares it
            if( !node.isMultipleAssignmentDeclaration() && PARAMS.equals(node.getVariableExpression().getName()) ) {
                node.getRightExpression().visit(this);
                shadowed = true;
                return;
            }
            super.visitDeclarationExpression(node);
        }

        @Override
        public void visitVariableExpression(VariableExpression node) {
            if( PARAMS.equals(node.getName()) && !shadowed )
                found = true;
        }

        private static boolean isShadowing(Parameter[] parameters) {
            return parameters != null && Arrays.stream(parameters).anyMatch(p -> PARAMS.equals(p.getName()));
        }
    }

    private static final String PIPELINE_NAME = "workflow";

    private static final String OUTPUT_NAME = "output";

    /**
     * The `params` and `output` blocks of an included pipeline can be included
     * as record types, so that a calling pipeline can declare a single param or
     * output for the entire pipeline instead of replicating each one.
     *
     * The params record type is *partial* -- every field is nullable, because
     * a param can be provided by the calling pipeline instead of the user, and
     * the pipeline validates its params when it is called.
     */
    private static ClassNode pipelineBlockType(ScriptNode sn, IncludeEntryNode entry) {
        if( "params".equals(entry.name) && sn.getParams() != null )
            return recordType(entry.getNameOrAlias(), List.of(sn.getParams().declarations), true);
        if( "output".equals(entry.name) && sn.getOutputs() != null )
            return recordType(entry.getNameOrAlias(), List.copyOf(sn.getOutputs().declarations), false);
        return null;
    }

    private static final ClassNode NULLABLE = ClassHelper.makeCached(Nullable.class);

    /**
     * Marks a record type synthesized for the `params` or `output` block of an
     * included pipeline, so that it can be distinguished from a definition of
     * the module with the same name.
     */
    public static final String PIPELINE_BLOCK_TYPE = "nextflow.pipelineBlockType";

    private static ClassNode recordType(String name, List<? extends Parameter> declarations, boolean nullable) {
        var cn = new RecordNode(name);
        cn.putNodeMetaData(PIPELINE_BLOCK_TYPE, Boolean.TRUE);
        for( var declaration : declarations ) {
            var fn = new FieldNode(declaration.getName(), java.lang.reflect.Modifier.PUBLIC, declaration.getType(), cn, null);
            fn.setDeclaringClass(cn);
            if( nullable )
                fn.addAnnotation(new AnnotationNode(NULLABLE));
            cn.addField(fn);
        }
        return cn;
    }

    private static void setPlaceholderTargets(IncludeNode node) {
        for( var entry : node.entries ) {
            if( entry.getTarget() == null ) {
                var target = new FunctionNode(entry.getNameOrAlias());
                entry.setTarget(target);
            }
        }
    }

    private URI getIncludeUri(URI base, String source) {
        if( ModuleResolver.isRemoteModule(source) ) {
            return RemoteModuleResolverProvider.getInstance()
                .resolve(source, projectDir)
                .normalize()
                .toUri();
        }
        else {
            var parent = Path.of(base).getParent();
            return getLocalIncludeUri(parent, source);
        }
    }

    private static URI getLocalIncludeUri(Path parent, String source) {
        Path includePath = parent.resolve(source);
        if( Files.isDirectory(includePath) )
            includePath = includePath.resolve("main.nf");
        else if( !source.endsWith(".nf") )
            includePath = Path.of(includePath.toString() + ".nf");
        return includePath.normalize().toUri();
    }

    private boolean isIncludeStale(IncludeNode node, URI includeUri) {
        if( changedUris == null || changedUris.contains(uri) || changedUris.contains(includeUri) )
            return true;
        for( var entry : node.entries ) {
            if( entry.getTarget() == null )
                return true;
        }
        return false;
    }

    private List<AnnotatedNode> getDefinitions(URI uri) {
        var scriptNode = (ScriptNode) compiler.getSource(uri).getAST();
        var result = new ArrayList<AnnotatedNode>();
        result.addAll(scriptNode.getWorkflows());
        result.addAll(scriptNode.getProcesses());
        result.addAll(scriptNode.getAgents());
        result.addAll(scriptNode.getFunctions());
        result.addAll(scriptNode.getTypes());
        return result;
    }

    /**
     * An entire pipeline -- the `params` / `workflow` / `output` trio of a
     * script -- can be included as a named workflow, using the `workflow`
     * keyword to refer to the entry workflow of the included script.
     */
    private static String definitionName(AnnotatedNode node) {
        if( node instanceof WorkflowNode wn && wn.isEntry() )
            return "workflow";
        return
            node instanceof ClassNode cn ? cn.getNameWithoutPackage() :
            node instanceof MethodNode mn ? mn.getName() :
            null;
    }

    @Override
    public void addError(String message, ASTNode node) {
        var cause = new ResolveIncludeError(message, node);
        var errorMessage = new SyntaxErrorMessage(cause, sourceUnit);
        errors.add(errorMessage);
    }

    public List<SyntaxErrorMessage> getErrors() {
        return errors;
    }

    public boolean isChanged() {
        return changed;
    }

    private class ResolveIncludeError extends SyntaxException implements PhaseAware {

        public ResolveIncludeError(String message, ASTNode node) {
            super(message, node);
        }

        @Override
        public int getPhase() {
            return Phases.INCLUDE_RESOLUTION;
        }
    }
}
