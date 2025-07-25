/*
 * Copyright 2013-2024, Seqera Labs
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

package nextflow.script.parser.v2;

import java.io.File;
import java.io.IOException;
import java.net.URI;
import java.nio.file.Path;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import groovy.lang.GroovyClassLoader;
import groovy.lang.GroovyCodeSource;
import com.google.common.hash.Hashing;
import nextflow.script.ast.WorkflowNode;
import nextflow.script.control.CallSiteCollector;
import nextflow.script.control.Compiler;
import nextflow.script.control.ModuleResolver;
import nextflow.script.control.OpCriteriaVisitor;
import nextflow.script.control.PathCompareVisitor;
import nextflow.script.control.ProcessNameResolver;
import nextflow.script.control.ResolveIncludeVisitor;
import nextflow.script.control.ScriptResolveVisitor;
import nextflow.script.control.ScriptToGroovyVisitor;
import nextflow.script.control.TypeCheckingVisitor;
import nextflow.script.parser.ScriptParserPluginFactory;
import org.codehaus.groovy.ast.ASTNode;
import org.codehaus.groovy.ast.ClassHelper;
import org.codehaus.groovy.ast.ClassNode;
import org.codehaus.groovy.ast.MethodNode;
import org.codehaus.groovy.control.CompilationFailedException;
import org.codehaus.groovy.control.CompilationUnit;
import org.codehaus.groovy.control.CompilerConfiguration;
import org.codehaus.groovy.control.Phases;
import org.codehaus.groovy.control.SourceUnit;
import org.codehaus.groovy.control.customizers.ImportCustomizer;
import org.codehaus.groovy.control.io.FileReaderSource;
import org.codehaus.groovy.control.io.ReaderSource;
import org.codehaus.groovy.control.io.StringReaderSource;
import org.codehaus.groovy.control.messages.SyntaxErrorMessage;

/**
 * Compile a Nextflow script into a Groovy class.
 *
 * @see groovy.lang.GroovyShell::parse()
 * @see groovy.lang.GroovyClassLoader::doParseClass()
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public class ScriptCompiler {

    private static final String DEFAULT_CODE_BASE = "/groovy/shell";
    private static final List<String> DEFAULT_IMPORT_NAMES = List.of(
        "java.nio.file.Path",
        "nextflow.Channel",
        "nextflow.util.Duration",
        "nextflow.util.MemoryUnit"
    );
    private static final String MAIN_CLASS_NAME = "Main";
    private static final String BASE_CLASS_NAME = "nextflow.script.BaseScript";

    private final CompilerConfiguration config;
    private final GroovyClassLoader loader;

    private Compiler compiler;

    public ScriptCompiler(boolean debug, Path targetDirectory, ClassLoader parent) {
        this(getConfig(debug, targetDirectory), parent);
    }

    public ScriptCompiler(CompilerConfiguration config, ClassLoader parent) {
        this.config = config;
        this.loader = new GroovyClassLoader(parent, config);
    }

    private static CompilerConfiguration getConfig(boolean debug, Path targetDirectory) {
        var importCustomizer = new ImportCustomizer();
        for ( var name : DEFAULT_IMPORT_NAMES )
            importCustomizer.addImports(name);
        importCustomizer.addImport("channel", "nextflow.Channel");
        importCustomizer.addStaticStars("nextflow.Nextflow");

        var config = new CompilerConfiguration();
        config.addCompilationCustomizers(importCustomizer);
        config.setScriptBaseClass(BASE_CLASS_NAME);
        config.setPluginFactory(new ScriptParserPluginFactory());
        config.setDebug(debug);
        if( targetDirectory != null )
            config.setTargetDirectory(targetDirectory.toFile());

        return config;
    }

    public CompileResult compile(String scriptText) throws IOException {
        return compile0(new GroovyCodeSource(scriptText, MAIN_CLASS_NAME, DEFAULT_CODE_BASE));
    }

    public CompileResult compile(File file) throws IOException {
        return compile0(new GroovyCodeSource(file, config.getSourceEncoding()));
    }

    public Collection<SourceUnit> getSources() {
        if( compiler == null )
            return null;
        return compiler.getSources().values();
    }

    public List<SyntaxErrorMessage> getErrors() {
        if( compiler == null )
            return null;
        return compiler.compilationUnit()
            .getErrorCollector()
            .getErrors()
            .stream()
            .map(e -> e instanceof SyntaxErrorMessage sem ? sem : null)
            .filter(sem -> sem != null)
            .toList();
    }

    private CompileResult compile0(GroovyCodeSource codeSource) throws IOException {
        // compile main script and included modules
        var unit = new ScriptCompilationUnit(config, loader);
        var su = codeSource.getFile() != null
            ? unit.createSourceUnit(MAIN_CLASS_NAME, codeSource.getFile())
            : unit.createSourceUnit(codeSource.getName(), codeSource.getScriptText());
        var collector = new ScriptClassLoader(loader).createCollector(unit, su);

        compiler = new Compiler(unit);
        compiler.addSource(su);

        unit.addSource(su);
        unit.setClassgenCallback(collector);
        unit.compile(Phases.CLASS_GENERATION);

        // collect script classes
        var classes = (List<Class>) collector.getLoadedClasses().stream()
            .map((o) -> 
                o instanceof Class c && BASE_CLASS_NAME.equals(c.getSuperclass().getName())
                    ? c
                    : null
            )
            .filter(c -> c != null)
            .toList();

        // extract main script class
        var main = classes.stream()
            .filter((c) -> MAIN_CLASS_NAME.equals(c.getSimpleName()))
            .findFirst()
            .get();

        var modules = collectModules(unit, classes);
        var processNames = new ProcessNameResolver(unit.getCallSites()).resolve(su);
        return new CompileResult(main, modules, processNames);
    }

    private Map<Path,Class> collectModules(ScriptCompilationUnit unit, List<Class> classes) {
        // match each module script class to the source path
        // using the class name
        var result = new HashMap<Path,Class>();
        for( var c : classes ) {
            for( var source : unit.getModules() ) {
                if( source.getName().equals(c.getSimpleName()) ) {
                    var path = Path.of(source.getSource().getURI());
                    result.put(path, c);
                    break;
                }
            }
        }
        return result;
    }

    public static record CompileResult(
        Class main,
        Map<Path,Class> modules,
        Set<String> processNames
    ) {}

    private static class ScriptClassLoader extends GroovyClassLoader {
        ScriptClassLoader(GroovyClassLoader parent) {
            super(parent);
        }

        @Override
        protected ClassCollector createCollector(CompilationUnit unit, SourceUnit su) {
            return super.createCollector(unit, su);
        }
    }

    private class ScriptCompilationUnit extends CompilationUnit {

        private static final List<ClassNode> DEFAULT_IMPORTS = defaultImports();

        private static List<ClassNode> defaultImports() {
            return DEFAULT_IMPORT_NAMES.stream()
                .map(name -> ClassHelper.makeWithoutCaching(name))
                .toList();
        }

        private SourceUnit entry;

        private Set<SourceUnit> modules;

        private Map<WorkflowNode, Map<String, MethodNode>> callSites = new IdentityHashMap<>();

        ScriptCompilationUnit(CompilerConfiguration configuration, GroovyClassLoader loader) {
            super(configuration, null, loader);
            super.addPhaseOperation(source -> analyze(source), Phases.CONVERSION);
        }

        Set<SourceUnit> getModules() {
            return modules;
        }

        Map<WorkflowNode, Map<String, MethodNode>> getCallSites() {
            return callSites;
        }

        @Override
        public void addPhaseOperation(ISourceUnitOperation op, int phase) {
            super.addPhaseOperation((source) -> {
                // skip main script on second conversion pass
                if( phase == Phases.CONVERSION && source == entry )
                    return;
                op.call(source);
            }, phase);
        }

        @Override
        public void addPhaseOperation(IPrimaryClassNodeOperation op, int phase) {
            super.addPhaseOperation((source, context, classNode) -> {
                // skip main script on second conversion pass
                if( phase == Phases.CONVERSION && source == entry )
                    return;
                op.call(source, context, classNode);
            }, phase);
        }

        private void analyze(SourceUnit source) {
            // on first pass, recursively add included modules to queue
            if( entry == null ) {
                modules = new ModuleResolver(compiler).resolve(source, uri -> createSourceUnit(uri));
                for( var su : modules )
                    addSource(su);
                entry = source;
                // wait for modules to be parsed before analyzing main script
                if( !modules.isEmpty() )
                    return;
            }

            // on second pass, all source files have been parsed
            // initialize script class
            var cn = source.getAST().getClasses().get(0);

            // perform strict syntax checking
            var includeResolver = new ResolveIncludeVisitor(source, compiler);
            includeResolver.visit();
            for( var error : includeResolver.getErrors() )
                source.getErrorCollector().addErrorAndContinue(error);
            new ScriptResolveVisitor(source, this, DEFAULT_IMPORTS, Collections.emptyList()).visit();
            if( source.getErrorCollector().hasErrors() )
                return;
            new TypeCheckingVisitor(source, false).visit();
            if( source.getErrorCollector().hasErrors() )
                return;

            // collect call sites for each workflow in the script
            callSites.putAll(new CallSiteCollector().apply(source));

            // convert to Groovy
            new ScriptToGroovyVisitor(source).visit();
            new PathCompareVisitor(source).visitClass(cn);
            new OpCriteriaVisitor(source).visitClass(cn);
        }

        SourceUnit createSourceUnit(URI uri) {
            return createSourceUnit(uniqueClassName(uri), new File(uri));
        }

        SourceUnit createSourceUnit(String name, File file) {
            return createSourceUnit(name, new FileReaderSource(file, getConfiguration()));
        }

        SourceUnit createSourceUnit(String name, String source) {
            return createSourceUnit(name, new StringReaderSource(source, getConfiguration()));
        }

        SourceUnit createSourceUnit(String name, ReaderSource source) {
            return new SourceUnit(
                    name,
                    source,
                    getConfiguration(),
                    getClassLoader(),
                    getErrorCollector());
        }

        /**
         * Create a unique name for a script class in order to avoid collisions
         * between scripts with the same base name.
         *
         * @param uri
         */
        private String uniqueClassName(URI uri) {
            var hash = Hashing
                    .sipHash24()
                    .newHasher()
                    .putUnencodedChars(uri.toString())
                    .hash();
            return "_nf_script_" + hash;
        }
    }

}
