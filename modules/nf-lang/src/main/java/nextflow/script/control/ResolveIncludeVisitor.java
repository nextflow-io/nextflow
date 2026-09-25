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
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import nextflow.script.ast.FunctionNode;
import nextflow.script.ast.IncludeEntryNode;
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
import org.codehaus.groovy.ast.FieldNode;
import org.codehaus.groovy.ast.MethodNode;
import org.codehaus.groovy.ast.Parameter;
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
            includeUri = ModuleResolver.getIncludeUri(uri, source, projectDir);
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
        var hasPipeline = false;
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
                hasPipeline |= PIPELINE_NAME.equals(includedName);
            }
            entry.setTarget(includedNode);
        }
        if( hasPipeline && !scriptNode.getParamsV1().isEmpty() )
            addError("An included pipeline cannot use legacy parameter declarations -- use the `params` block instead", node);
    }

    private static boolean isEntryWorkflow(String name, AnnotatedNode node) {
        return PIPELINE_NAME.equals(name) && node instanceof WorkflowNode wn && wn.isEntry();
    }

    private static final String PIPELINE_NAME = "workflow";

    /**
     * The `params` and `output` blocks of an included pipeline can be included
     * as record types, so that a calling pipeline can refer to the params or
     * outputs of the pipeline as a whole instead of replicating each one.
     *
     * The params record type is *partial* -- every field is nullable, because
     * a param can be provided by the calling pipeline instead of the user, and
     * the pipeline validates its params when it is called.
     */
    private static ClassNode pipelineBlockType(ScriptNode sn, IncludeEntryNode entry) {
        if( "params".equals(entry.name) && sn.getParams() != null )
            return recordType(entry.getNameOrAlias(), sn.getParams(), List.of(sn.getParams().declarations), true);
        if( "output".equals(entry.name) && sn.getOutputs() != null )
            return recordType(entry.getNameOrAlias(), sn.getOutputs(), List.copyOf(sn.getOutputs().declarations), false);
        return null;
    }

    private static final ClassNode NULLABLE = ClassHelper.makeCached(Nullable.class);

    /**
     * Marks a record type synthesized for the `params` or `output` block of an
     * included pipeline, so that it can be distinguished from a definition of
     * the module with the same name. The metadata value is the block that the
     * type was synthesized from, so that tooling can navigate to it.
     */
    public static final String PIPELINE_BLOCK_TYPE = "nextflow.pipelineBlockType";

    /**
     * Get the `params` or `output` block that a record type was synthesized
     * from, or null if the type is a definition of the module.
     *
     * @param cn
     */
    public static ASTNode getPipelineBlock(ClassNode cn) {
        return (ASTNode) cn.getNodeMetaData(PIPELINE_BLOCK_TYPE);
    }

    private static ClassNode recordType(String name, ASTNode block, List<? extends Parameter> declarations, boolean nullable) {
        var cn = new RecordNode(name);
        cn.putNodeMetaData(PIPELINE_BLOCK_TYPE, block);
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
