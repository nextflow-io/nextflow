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

import java.io.File;
import java.net.URI;
import java.nio.file.Path;
import java.security.CodeSource;
import java.util.HashMap;
import java.util.Map;

import groovy.lang.GroovyClassLoader;
import nextflow.script.ast.IncludeNode;
import org.antlr.v4.runtime.RecognitionException;
import org.codehaus.groovy.control.CompilationFailedException;
import org.codehaus.groovy.control.CompilationUnit;
import org.codehaus.groovy.control.CompilerConfiguration;
import org.codehaus.groovy.control.ErrorCollector;
import org.codehaus.groovy.control.SourceUnit;

/**
 * Compiler that can lookup source units by URI.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public class Compiler {

    private CompilationUnit compilationUnit;

    private Map<URI, SourceUnit> sourcesByUri = new HashMap<>();

    private Map<String, RuntimeException> includeFailures = new HashMap<>();

    public Compiler(CompilerConfiguration configuration, GroovyClassLoader classLoader) {
        this(new CompilationUnit(configuration, null, classLoader));
    }

    public Compiler(CompilationUnit compilationUnit) {
        this.compilationUnit = compilationUnit;
    }

    public CompilationUnit compilationUnit() {
        return compilationUnit;
    }

    protected CompilerConfiguration configuration() {
        return compilationUnit().getConfiguration();
    }

    protected GroovyClassLoader classLoader() {
        return compilationUnit().getClassLoader();
    }

    public SourceUnit createSourceUnit(File file) {
        return new SourceUnit(
            file,
            configuration(),
            classLoader(),
            createErrorCollector());
    }

    public SourceUnit createSourceUnit(String name, String contents) {
        return new SourceUnit(
            name,
            contents,
            configuration(),
            classLoader(),
            createErrorCollector());
    }

    protected ErrorCollector createErrorCollector() {
        return new LazyErrorCollector(configuration());
    }

    public void addSource(SourceUnit source) {
        sourcesByUri.put(source.getSource().getURI(), source);
    }

    public Map<URI, SourceUnit> getSources() {
        return sourcesByUri;
    }

    public SourceUnit getSource(URI uri) {
        return sourcesByUri.get(uri);
    }

    /**
     * Resolve an include URI, remembering failures so that a failing remote include
     * (e.g. a pipeline auto-install) is attempted once per compiler.
     */
    public URI getIncludeUri(URI uri, IncludeNode node, Path projectDir) {
        var key = uri + "|" + node.source.getText();
        var failure = includeFailures.get(key);
        if( failure != null )
            throw failure;
        try {
            return ModuleResolver.getIncludeUri(uri, node, projectDir);
        }
        catch( RuntimeException e ) {
            includeFailures.put(key, e);
            throw e;
        }
    }

    public void compile(SourceUnit source) {
        try {
            source.parse();
            source.buildAST();
        }
        catch( RecognitionException e ) {
        }
        catch( CompilationFailedException e ) {
        }
    }

}
