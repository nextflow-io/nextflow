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

package nextflow.script.parser.v2

import java.nio.file.Path

import com.google.common.hash.Hashing
import groovy.transform.CompileStatic
import nextflow.Session
import nextflow.ast.NextflowXformImpl
import nextflow.ast.OpXformImpl
import nextflow.script.BaseScript
import nextflow.script.control.Compiler
import nextflow.script.control.ModuleResolver
import nextflow.script.control.ResolveIncludeVisitor
import nextflow.script.control.ScriptResolveVisitor
import nextflow.script.control.TypeCheckingVisitor
import nextflow.script.parser.ScriptParserPluginFactory
import org.codehaus.groovy.ast.ASTNode
import org.codehaus.groovy.ast.ClassNode
import org.codehaus.groovy.control.CompilationFailedException
import org.codehaus.groovy.control.CompilationUnit
import org.codehaus.groovy.control.CompilerConfiguration
import org.codehaus.groovy.control.Phases
import org.codehaus.groovy.control.SourceUnit
import org.codehaus.groovy.control.customizers.ImportCustomizer
import org.codehaus.groovy.control.io.FileReaderSource
import org.codehaus.groovy.control.io.ReaderSource
import org.codehaus.groovy.control.io.StringReaderSource
import org.codehaus.groovy.control.messages.SyntaxErrorMessage

/**
 * Compile a Nextflow script into a Groovy class.
 *
 * @see groovy.lang.GroovyShell::parse()
 * @see groovy.lang.GroovyClassLoader::doParseClass()
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
public class ScriptCompiler {

    public static final String DEFAULT_CODE_BASE = "/groovy/shell"
    public static final String MAIN_CLASS_NAME = "Main"

    private final CompilerConfiguration config
    private final GroovyClassLoader loader

    private Compiler compiler

    public ScriptCompiler(Session session) {
        this(getConfig(session), session.classLoader)
    }

    public ScriptCompiler(CompilerConfiguration config, ClassLoader parent) {
        this.config = config
        this.loader = new GroovyClassLoader(parent, config)
    }

    private static CompilerConfiguration getConfig(Session session) {
        final importCustomizer = new ImportCustomizer()
        importCustomizer.addImports( java.nio.file.Path.name )
        importCustomizer.addImports( nextflow.Channel.name )
        importCustomizer.addImports( nextflow.util.Duration.name )
        importCustomizer.addImports( nextflow.util.MemoryUnit.name )
        importCustomizer.addImport( "channel", nextflow.Channel.name )
        importCustomizer.addStaticStars( nextflow.Nextflow.name )

        final config = new CompilerConfiguration()
        config.addCompilationCustomizers( importCustomizer )
        config.setScriptBaseClass(BaseScript.class.getName())
        config.setPluginFactory(new ScriptParserPluginFactory())

        if( session?.debug )
            config.debug = true

        if( session?.classesDir )
            config.setTargetDirectory(session.classesDir.toFile())

        return config
    }

    public CompileResult compile(String scriptText) {
        return compile0(new GroovyCodeSource(scriptText, MAIN_CLASS_NAME, DEFAULT_CODE_BASE))
    }

    public CompileResult compile(File file) {
        return compile0(new GroovyCodeSource(file, config.getSourceEncoding()))
    }

    public Collection<SourceUnit> getSources() {
        if( !compiler )
            return null
        return compiler.getSources().values()
    }

    public List<SyntaxErrorMessage> getErrors() {
        if( !compiler )
            return null
        return compiler.compilationUnit()
            .getErrorCollector()
            .getErrors()
            .stream()
            .filter(e -> e instanceof SyntaxErrorMessage)
            .map(e -> (SyntaxErrorMessage) e)
            .toList()
    }

    private CompileResult compile0(GroovyCodeSource codeSource) {
        // compile main script and included modules
        final unit = new ScriptCompilationUnit(config, loader)
        final su = codeSource.getFile()
            ? unit.getSource(MAIN_CLASS_NAME, codeSource.getFile())
            : unit.getSource(codeSource.getName(), codeSource.getScriptText())
        final collector = new ScriptClassLoader(loader).createCollector(unit, su)

        compiler = new Compiler(unit)
        compiler.addSource(su)

        unit.addSource(su)
        unit.setClassgenCallback(collector)
        unit.compile(Phases.CLASS_GENERATION)

        // collect script classes
        final classes = collector.getLoadedClasses().stream()
            .map((o) -> (Class) o)
            .filter((c) -> c.getSuperclass() == BaseScript.class)
            .toList()

        // extract main script class
        final main = classes.stream()
            .filter(c -> c.getSimpleName() == MAIN_CLASS_NAME)
            .findFirst()
            .get()

        // match each module script class to the source path
        // using the class name
        final modules = new HashMap<Path,Class>()
        for( final c : classes ) {
            for( final source : unit.getModules() ) {
                if( source.getName() == c.getSimpleName() ) {
                    final path = Path.of(source.getSource().getURI())
                    modules.put(path, c)
                    break
                }
            }
        }
        return new CompileResult(main, modules)
    }

    public static record CompileResult(
        Class main,
        Map<Path,Class> modules
    ) {}

    private static class ScriptClassLoader extends GroovyClassLoader {
        public ScriptClassLoader(GroovyClassLoader parent) {
            super(parent)
        }

        @Override
        public ClassCollector createCollector(CompilationUnit unit, SourceUnit su) {
            return super.createCollector(unit, su)
        }
    }

    private class ScriptCompilationUnit extends CompilationUnit {

        private static final List<ClassNode> DEFAULT_IMPORTS = List.of(
            new ClassNode( java.nio.file.Path ),
            new ClassNode( nextflow.Channel ),
            new ClassNode( nextflow.util.Duration ),
            new ClassNode( nextflow.util.MemoryUnit ),
        )

        private SourceUnit entry

        private Set<SourceUnit> modules

        public ScriptCompilationUnit(CompilerConfiguration configuration, GroovyClassLoader loader) {
            super(configuration, null, loader)
            super.addPhaseOperation(source -> analyze(source), Phases.CONVERSION)
        }

        public Set<SourceUnit> getModules() {
            return modules
        }

        @Override
        public void addPhaseOperation(final ISourceUnitOperation op, final int phase) {
            super.addPhaseOperation((source) -> {
                // skip main script on second conversion pass
                if( phase == Phases.CONVERSION && source == entry )
                    return
                op.call(source)
            }, phase)
        }

        @Override
        public void addPhaseOperation(final IPrimaryClassNodeOperation op, final int phase) {
            super.addPhaseOperation((source, context, classNode) -> {
                // skip main script on second conversion pass
                if( phase == Phases.CONVERSION && source == entry )
                    return
                op.call(source, context, classNode)
            }, phase)
        }

        private void analyze(SourceUnit source) {
            // on first pass, recursively add included modules to queue
            if( entry == null ) {
                modules = new ModuleResolver(compiler).resolve(source, uri -> getSource(uri))
                for( final su : modules )
                    addSource(su)
                entry = source
                // wait for modules to be parsed before analyzing main script
                if( !modules.isEmpty() )
                    return
            }

            // on second pass, all source files have been parsed
            // initialize script class
            final cn = source.getAST().getClasses().get(0)

            // perform strict syntax checking
            final includeResolver = new ResolveIncludeVisitor(source, compiler, Collections.emptySet())
            includeResolver.visit()
            for( final error : includeResolver.getErrors() )
                source.getErrorCollector().addErrorAndContinue(error)
            new ScriptResolveVisitor(source, this, DEFAULT_IMPORTS, Collections.emptyList()).visit()
            if( source.getErrorCollector().hasErrors() )
                return
            new TypeCheckingVisitor(source, false).visit()
            if( source.getErrorCollector().hasErrors() )
                return

            // convert to Groovy
            final astNodes = new ASTNode[] { cn, cn }
            new ScriptToGroovyVisitor(source).visit()
            new NextflowXformImpl().visit(astNodes, source)
            new OpXformImpl().visit(astNodes, source)
        }

        SourceUnit getSource(URI uri) {
            return getSource(uniqueClassName(uri), new File(uri))
        }

        SourceUnit getSource(String name, File file) {
            return getSource(name, new FileReaderSource(file, getConfiguration()))
        }

        SourceUnit getSource(String name, String source) {
            return getSource(name, new StringReaderSource(source, getConfiguration()))
        }

        SourceUnit getSource(String name, ReaderSource source) {
            return new SourceUnit(
                    name,
                    source,
                    getConfiguration(),
                    getClassLoader(),
                    getErrorCollector())
        }

        /**
         * Create a unique name for a script class in order to avoid collisions
         * between scripts with the same base name.
         *
         * @param uri
         */
        private String uniqueClassName(URI uri) {
            final hash = Hashing
                    .sipHash24()
                    .newHasher()
                    .putUnencodedChars(uri.toString())
                    .hash()
            return "_nf_script_$hash"
        }
    }

}
