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
import org.codehaus.groovy.ast.ASTNode;
import org.codehaus.groovy.ast.AnnotatedNode;
import org.codehaus.groovy.ast.ClassNode;
import org.codehaus.groovy.ast.FieldNode;
import org.codehaus.groovy.ast.MethodNode;
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
            // a `params` entry that doesn't match a definition of the module
            // refers to the params block of the pipeline
            var target = definitions.stream()
                .filter(defNode -> includedName.equals(definitionName(defNode)))
                .findFirst()
                .orElseGet(() -> paramsBlockType(scriptNode, entry));
            if( target == null ) {
                addError("Included name '" + includedName + "' is not defined in module '" + includeUri.getPath() + "'", node);
                continue;
            }
            hasPipeline |= target instanceof WorkflowNode wn && wn.isEntry()
                || target instanceof ClassNode cn && ScriptNode.isPipelineParams(cn);
            entry.setTarget(target);
        }
        if( hasPipeline && !((ScriptNode) sourceUnit.getAST()).isTypingEnabled() )
            addError("Including a pipeline requires `nextflow.enable.types = true` in the including script", node);
        if( hasPipeline && !scriptNode.isTypingEnabled() )
            addError("An included pipeline must enable static typing -- set `nextflow.enable.types = true` in '" + includeUri.getPath() + "'", node);
    }

    /**
     * Synthesize a partial record type (all fields nullable) from
     * the params block of an included pipeline.
     */
    private static ClassNode paramsBlockType(ScriptNode sn, IncludeEntryNode entry) {
        var block = sn.getParams();
        if( !"params".equals(entry.name) || block == null )
            return null;
        var cn = new RecordNode(entry.getNameOrAlias());
        ScriptNode.setPipelineParams(cn);
        for( var declaration : block.declarations ) {
            cn.addField(new FieldNode(declaration.getName(), Modifier.PUBLIC, declaration.getType(), cn, null));
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
