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

import java.security.CodeSource

import groovy.transform.CompileStatic
import nextflow.Session
import nextflow.ast.NextflowXformImpl
import nextflow.ast.OpXformImpl
import nextflow.script.BaseScript
import nextflow.script.control.ResolveVisitor
import nextflow.script.parser.ScriptParserPluginFactory
import org.codehaus.groovy.ast.ASTNode
import org.codehaus.groovy.ast.ClassNode
import org.codehaus.groovy.control.CompilationFailedException
import org.codehaus.groovy.control.CompilationUnit
import org.codehaus.groovy.control.CompilerConfiguration
import org.codehaus.groovy.control.Phases
import org.codehaus.groovy.control.customizers.ImportCustomizer
import org.codehaus.groovy.runtime.InvokerHelper

/**
 * Compile a Nextflow script into a Groovy class.
 *
 * @see groovy.lang.GroovyShell
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
public class ScriptCompiler {

    public static final String DEFAULT_CODE_BASE = "/groovy/shell"

    private final Binding context
    private final CompilerConfiguration config
    private final ScriptClassLoader loader

    public ScriptCompiler(Session session, Binding binding) {
        this(session.classLoader, binding, getConfig(session))
    }

    public ScriptCompiler(ClassLoader parent, Binding binding, CompilerConfiguration config) {
        this.loader = new ScriptClassLoader(parent, config)
        this.context = binding
        this.config = config
    }

    public Script compile(String scriptText, String fileName) throws CompilationFailedException {
        GroovyCodeSource gcs = new GroovyCodeSource(scriptText, fileName, DEFAULT_CODE_BASE)
        return InvokerHelper.createScript(parseClass(gcs), context)
    }

    public Script compile(File file) throws CompilationFailedException, IOException {
        GroovyCodeSource gcs = new GroovyCodeSource(file, config.getSourceEncoding())
        return InvokerHelper.createScript(parseClass(gcs), context)
    }

    private Class parseClass(GroovyCodeSource codeSource) throws CompilationFailedException {
        // Don't cache scripts
        return loader.parseClass(codeSource, false)
    }

    private static CompilerConfiguration getConfig(Session session) {
        final importCustomizer = new ImportCustomizer()
        importCustomizer.addImports( java.nio.file.Path.name )
        importCustomizer.addImports( nextflow.Channel.name )
        importCustomizer.addImports( nextflow.util.Duration.name )
        importCustomizer.addImports( nextflow.util.MemoryUnit.name )
        importCustomizer.addImport( 'channel', nextflow.Channel.name )
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

    private static class ScriptClassLoader extends GroovyClassLoader {

        public ScriptClassLoader(ClassLoader loader, CompilerConfiguration config) {
            super(loader, config)
        }

        @Override
        protected CompilationUnit createCompilationUnit(CompilerConfiguration config, CodeSource source) {
            return new ScriptCompilationUnit(config, source, this)
        }
    }

    private static class ScriptCompilationUnit extends CompilationUnit {

        private static final List<ClassNode> DEFAULT_IMPORTS = List.of(
            new ClassNode( java.nio.file.Path ),
            new ClassNode( nextflow.Channel ),
            new ClassNode( nextflow.util.Duration ),
            new ClassNode( nextflow.util.MemoryUnit ),
        )

        public ScriptCompilationUnit(CompilerConfiguration configuration, CodeSource codeSource, GroovyClassLoader loader) {
            super(configuration, codeSource, loader)

            addPhaseOperation(source -> {
                // initialize script class
                final cn = source.getAST().getClasses().get(0)
                final astNodes = new ASTNode[] { cn, cn }

                new ResolveVisitor(source, this, DEFAULT_IMPORTS, new ArrayList<>()).visit()
                new ScriptToGroovyVisitor(source).visit()
                new NextflowXformImpl().visit(astNodes, source)
                new OpXformImpl().visit(astNodes, source)
            }, Phases.CONVERSION)
        }
    }

}
