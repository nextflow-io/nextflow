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
import java.util.LinkedList;
import java.util.HashSet;
import java.util.Set;
import java.util.function.Function;

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
        var uri = sourceUnit.getSource().getURI();
        var includeUri = getIncludeUri(uri, source);
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

    private static URI getIncludeUri(URI uri, String source) {
        Path includePath = Path.of(uri).getParent().resolve(source);
        if( Files.isDirectory(includePath) )
            includePath = includePath.resolve("main.nf");
        else if( !source.endsWith(".nf") )
            includePath = Path.of(includePath.toString() + ".nf");
        return includePath.normalize().toUri();
    }

}
