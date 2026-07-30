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
import java.util.List;
import java.util.Set;

import nextflow.module.spi.RemoteModuleResolverProvider;
import nextflow.script.ast.FunctionNode;
import nextflow.script.ast.IncludeNode;
import nextflow.script.ast.ScriptNode;
import nextflow.script.ast.ScriptVisitorSupport;
import org.codehaus.groovy.ast.ASTNode;
import org.codehaus.groovy.ast.AnnotatedNode;
import org.codehaus.groovy.ast.ClassNode;
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
            includeUri = getIncludeUri(source);
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
        var definitions = getDefinitions(includeUri);
        for( var entry : node.entries ) {
            var includedName = entry.name;
            var includedNode = definitions.stream()
                .filter(defNode -> includedName.equals(definitionName(defNode)))
                .findFirst();
            if( !includedNode.isPresent() ) {
                addError("Included name '" + includedName + "' is not defined in module '" + includeUri.getPath() + "'", node);
                continue;
            }
            entry.setTarget(includedNode.get());
        }
    }

    private static void setPlaceholderTargets(IncludeNode node) {
        for( var entry : node.entries ) {
            if( entry.getTarget() == null ) {
                var target = new FunctionNode(entry.getNameOrAlias());
                entry.setTarget(target);
            }
        }
    }

    private URI getIncludeUri(String source) {
        if( ModuleResolver.isRemoteModule(source) ) {
            return RemoteModuleResolverProvider.getInstance()
                .resolve(source, projectDir)
                .normalize()
                .toUri();
        }
        else {
            var parent = Path.of(uri).getParent();
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
        result.addAll(scriptNode.getFunctions());
        result.addAll(scriptNode.getTypes());
        return result;
    }

    private static String definitionName(AnnotatedNode node) {
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
