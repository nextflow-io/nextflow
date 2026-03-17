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
import java.util.LinkedList;
import java.util.HashSet;
import java.util.Set;
import java.util.function.Function;

import nextflow.module.spi.RemoteModuleResolverProvider;
import nextflow.script.ast.IncludeNode;
import nextflow.script.ast.ScriptNode;
import org.codehaus.groovy.control.SourceUnit;

/**
 * Resolve and compile all modules included (directly or indirectly)
 * by the main script.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public class ModuleResolver {

    private Compiler compiler;

    public ModuleResolver(Compiler compiler) {
        this.compiler = compiler;
    }

    /**
     * Resolve all modules included by a script.
     *
     * @param entry the main script
     * @param sourceResolver function that generates the source unit for a given file
     */
    public Set<SourceUnit> resolve(SourceUnit entry, Function<URI,SourceUnit> sourceResolver) {
        var modules = new HashSet<SourceUnit>();
        var queuedSources = new LinkedList<SourceUnit>();
        compiler.addSource(entry);
        queuedSources.add(entry);
        while( !queuedSources.isEmpty() ) {
            var source = queuedSources.remove();
            if( source.getAST() == null )
                continue;
            var sn = (ScriptNode) source.getAST();
            for( var in : sn.getIncludes() ) {
                var includeSource = resolveInclude(in, source, sourceResolver);
                if( includeSource == null )
                    continue;
                modules.add(includeSource);
                queuedSources.add(includeSource);
            }
        }
        return modules;
    }

    private SourceUnit resolveInclude(IncludeNode node, SourceUnit sourceUnit, Function<URI,SourceUnit> sourceResolver) {
        var source = node.source.getText();
        if( source.startsWith("plugin/") )
            return null;

        var includeUri = isRemoteModule(source) ?
            RemoteModuleResolverProvider.getInstance().resolve(source).normalize().toUri() :
            getIncludeUri(Path.of(sourceUnit.getSource().getURI()).getParent(), source);
        if( compiler.getSource(includeUri) != null )
            return null;
        if( !Files.exists(Path.of(includeUri)) )
            return null;
        var includeSource = sourceResolver.apply(includeUri);
        compiler.addSource(includeSource);
        compiler.compile(includeSource);
        if( includeSource.getAST() == null )
            return null;
        return includeSource;
    }

    /**
     * Module name pattern matching the canonical format used by ModuleReference.
     * Scope: lowercase alphanumeric with dots/underscores/hyphens.
     * Name: one or more slash-separated segments, each lowercase alphanumeric with dots/underscores/hyphens.
     */
    static final String REMOTE_MODULE_PATTERN = "^[a-z0-9][a-z0-9._\\-]*/[a-z][a-z0-9._\\-]*(/[a-z][a-z0-9._\\-]*)*$";

    static boolean isRemoteModule(String source) {
        if( source.startsWith("/") || source.startsWith("./") || source.startsWith("../") )
            return false;
        return source.matches(REMOTE_MODULE_PATTERN);
    }

    private static URI getIncludeUri(Path parent, String source) {
        Path includePath = parent.resolve(source);
        if( Files.isDirectory(includePath) )
            includePath = includePath.resolve("main.nf");
        else if( !source.endsWith(".nf") )
            includePath = Path.of(includePath.toString() + ".nf");
        return includePath.normalize().toUri();
    }

}
