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

import java.net.URI;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import nextflow.script.ast.FunctionNode;
import nextflow.script.ast.IncludeNode;
import nextflow.script.ast.ScriptNode;
import nextflow.script.ast.ScriptVisitorSupport;
import org.codehaus.groovy.ast.ASTNode;
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

    private Compiler compiler;

    private Set<URI> changedUris;

    private List<SyntaxErrorMessage> errors = new ArrayList<>();

    private boolean changed;

    public ResolveIncludeVisitor(SourceUnit sourceUnit, Compiler compiler, Set<URI> changedUris) {
        this.sourceUnit = sourceUnit;
        this.uri = sourceUnit.getSource().getURI();
        this.compiler = compiler;
        this.changedUris = changedUris;
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
        var includeUri = getIncludeUri(uri, source);
        if( !isIncludeStale(node, includeUri) )
            return;
        changed = true;
        for( var module : node.modules )
            module.setTarget(null);
        var includeUnit = compiler.getSource(includeUri);
        if( includeUnit == null ) {
            addError("Invalid include source: '" + includeUri + "'", node);
            return;
        }
        if( includeUnit.getAST() == null ) {
            addError("Module could not be parsed: '" + includeUri + "'", node);
            return;
        }
        var definitions = getDefinitions(includeUri);
        for( var module : node.modules ) {
            var includedName = module.name;
            var includedNode = definitions.stream()
                .filter(defNode -> includedName.equals(defNode.getName()))
                .findFirst();
            if( !includedNode.isPresent() ) {
                addError("Included name '" + includedName + "' is not defined in module '" + includeUri + "'", node);
                continue;
            }
            module.setTarget(includedNode.get());
        }
    }

    private static void setPlaceholderTargets(IncludeNode node) {
        for( var module : node.modules ) {
            if( module.getTarget() == null ) {
                var target = new FunctionNode(module.getNameOrAlias());
                module.setTarget(target);
            }
        }
    }

    private static URI getIncludeUri(URI uri, String source) {
        Path includePath = Path.of(uri).getParent().resolve(source);
        if( Files.isDirectory(includePath) )
            includePath = includePath.resolve("main.nf");
        else if( !source.endsWith(".nf") )
            includePath = Path.of(includePath.toString() + ".nf");
        return includePath.normalize().toUri();
    }

    private boolean isIncludeStale(IncludeNode node, URI includeUri) {
        if( changedUris.contains(uri) || changedUris.contains(includeUri) )
            return true;
        for( var module : node.modules ) {
            if( module.getTarget() == null )
                return true;
        }
        return false;
    }

    private List<MethodNode> getDefinitions(URI uri) {
        var scriptNode = (ScriptNode) compiler.getSource(uri).getAST();
        var result = new ArrayList<MethodNode>();
        result.addAll(scriptNode.getWorkflows());
        result.addAll(scriptNode.getProcesses());
        result.addAll(scriptNode.getFunctions());
        return result;
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
