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
import java.util.Collection;
import java.util.LinkedList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Function;

import nextflow.module.spi.RemoteModuleResolver;
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
    private Path projectDir;

    public ModuleResolver(Path projectDir, Compiler compiler) {
        this.compiler = compiler;
        this.projectDir = projectDir;
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

    /**
     * Order the given sources so that a script appears after every module
     * it includes, via depth-first post-order traversal. Only includes
     * among the given sources are considered.
     *
     * @param sources
     */
    public List<SourceUnit> orderByDependencies(Collection<SourceUnit> sources) {
        var byUri = new HashMap<URI,SourceUnit>();
        for( var source : sources )
            byUri.put(source.getSource().getURI(), source);
        var ordered = new ArrayList<SourceUnit>(sources.size());
        var visited = new HashSet<URI>();
        for( var source : sources )
            visitDependencies(source, byUri, visited, ordered);
        return ordered;
    }

    private void visitDependencies(SourceUnit source, Map<URI,SourceUnit> byUri, Set<URI> visited, List<SourceUnit> ordered) {
        var uri = source.getSource().getURI();
        if( !visited.add(uri) )
            return;
        if( source.getAST() instanceof ScriptNode sn ) {
            for( var in : sn.getIncludes() ) {
                var include = in.source.getText();
                if( include.startsWith("plugin/") )
                    continue;
                var dep = byUri.get(getIncludeUri(uri, include, projectDir));
                if( dep != null )
                    visitDependencies(dep, byUri, visited, ordered);
            }
        }
        ordered.add(source);
    }

    private SourceUnit resolveInclude(IncludeNode node, SourceUnit sourceUnit, Function<URI,SourceUnit> sourceResolver) {
        var source = node.source.getText();
        if( source.startsWith("plugin/") )
            return null;

        var uri = sourceUnit.getSource().getURI();
        var includeUri = getIncludeUri(uri, source, projectDir);
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
     * @return true if the given include source refers to a local module, i.e. it is a path to a
     * script. Any other include source is a remote module reference -- a malformed one is
     * reported as an invalid module reference by the resolver.
     */
    public static boolean isLocalModule(String source) {
        return source.startsWith("/") || source.startsWith("./") || source.startsWith("../");
    }

    /**
     * Resolve an include source to the URI of the included script.
     *
     * @param uri the URI of the including script
     * @param source the include source
     * @param projectDir the project directory, used to resolve remote modules
     */
    public static URI getIncludeUri(URI uri, String source, Path projectDir) {
        var localUri = getLocalIncludeUri(uri, source);
        if( localUri != null )
            return localUri;
        // Resolve a remote module relative to the including module's directory
        // (context-relative), so a workflow module's own dependencies are found under its
        // nested `modules/` directory (nested vendoring). Any other script -- the entry
        // script, or a plain local script -- resolves against the project directory.
        var base = RemoteModuleResolver.resolveBaseDir(uri, projectDir);
        return RemoteModuleResolverProvider.getInstance()
            .resolve(source, base)
            .normalize()
            .toUri();
    }

    /**
     * Resolve a local include source to the URI of the included script.
     *
     * @param uri the URI of the including script
     * @param source the include source
     * @return the include URI, or null if the source is not a local module
     */
    public static URI getLocalIncludeUri(URI uri, String source) {
        if( !isLocalModule(source) )
            return null;
        Path includePath = Path.of(uri).getParent().resolve(source);
        if( Files.isDirectory(includePath) )
            includePath = includePath.resolve("main.nf");
        else if( !source.endsWith(".nf") )
            includePath = Path.of(includePath.toString() + ".nf");
        return includePath.normalize().toUri();
    }

}
