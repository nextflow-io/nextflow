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

package nextflow.config.parser.v2;

import java.io.IOException;
import java.nio.file.Path;
import java.util.List;

import com.google.common.hash.Hashing;
import groovy.lang.Binding;
import groovy.lang.GroovyClassLoader;
import groovy.lang.Script;
import nextflow.config.control.ConfigResolveVisitor;
import nextflow.config.control.ConfigToGroovyVisitor;
import nextflow.config.control.ResolveIncludeVisitor;
import nextflow.config.control.StringReaderSourceWithURI;
import nextflow.config.control.StripSecretsVisitor;
import nextflow.config.parser.ConfigParserPluginFactory;
import nextflow.script.control.Compiler;
import nextflow.script.control.PathCompareVisitor;
import org.codehaus.groovy.ast.ASTNode;
import org.codehaus.groovy.ast.ClassHelper;
import org.codehaus.groovy.ast.ClassNode;
import org.codehaus.groovy.control.CompilationUnit;
import org.codehaus.groovy.control.CompilerConfiguration;
import org.codehaus.groovy.control.Phases;
import org.codehaus.groovy.control.SourceUnit;
import org.codehaus.groovy.control.customizers.ImportCustomizer;
import org.codehaus.groovy.control.io.StringReaderSource;
import org.codehaus.groovy.control.messages.SyntaxErrorMessage;
import org.codehaus.groovy.runtime.InvokerHelper;

/**
 * Compile a Nextflow config file into a Groovy class.
 *
 * @see groovy.lang.GroovyShell::parse()
 * @see groovy.lang.GroovyClassLoader::doParseClass()
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public class ConfigCompiler {

    private static final List<String> DEFAULT_IMPORT_NAMES = List.of(
        "nextflow.util.Duration",
        "nextflow.util.MemoryUnit"
    );
    private static final String BASE_CLASS_NAME = "nextflow.config.parser.v2.ConfigDsl";

    private final CompilerConfiguration config;

    private final GroovyClassLoader loader;

    private boolean renderClosureAsString;

    private boolean stripSecrets;

    private Compiler compiler;

    private SourceUnit sourceUnit;

    public ConfigCompiler(boolean renderClosureAsString, boolean stripSecrets) {
        this.config = getConfig();
        this.loader = new GroovyClassLoader(new GroovyClassLoader(), config);
        this.renderClosureAsString = renderClosureAsString;
        this.stripSecrets = stripSecrets;
    }

    private static CompilerConfiguration getConfig() {
        var importCustomizer = new ImportCustomizer();
        for ( var name : DEFAULT_IMPORT_NAMES )
            importCustomizer.addImports(name);

        var config = new CompilerConfiguration();
        config.addCompilationCustomizers(importCustomizer);
        config.setScriptBaseClass(BASE_CLASS_NAME);
        config.setPluginFactory(new ConfigParserPluginFactory());

        return config;
    }

    public Script compile(String text, Path path) throws IOException {
        return compile0(text, path);
    }

    private Script compile0(String text, Path path) throws IOException {
        // prepare compiler and source unit
        var unit = new ConfigCompilationUnit(config, loader);
        sourceUnit = unit.createSourceUnit(text, path);
        var collector = new ConfigClassLoader(loader).createCollector(unit, sourceUnit);

        compiler = new Compiler(unit);
        compiler.addSource(sourceUnit);

        // compile config file
        unit.addSource(sourceUnit);
        unit.setClassgenCallback(collector);
        unit.compile(Phases.CLASS_GENERATION);

        // return compiled script class
        var clazz = (Class) collector.getLoadedClasses().stream().findFirst().orElse(null);
        return InvokerHelper.createScript(clazz, new Binding());
    }

    public SourceUnit getSource() {
        return sourceUnit;
    }

    public List<SyntaxErrorMessage> getErrors() {
        if( sourceUnit == null )
            return null;
        return sourceUnit
            .getErrorCollector()
            .getErrors()
            .stream()
            .map(e -> e instanceof SyntaxErrorMessage sem ? sem : null)
            .filter(sem -> sem != null)
            .toList();
    }

    private static class ConfigClassLoader extends GroovyClassLoader {
        ConfigClassLoader(GroovyClassLoader parent) {
            super(parent);
        }

        @Override
        protected ClassCollector createCollector(CompilationUnit unit, SourceUnit su) {
            return super.createCollector(unit, su);
        }
    }

    private class ConfigCompilationUnit extends CompilationUnit {

        private static final List<ClassNode> DEFAULT_IMPORTS = defaultImports();

        private static List<ClassNode> defaultImports() {
            return DEFAULT_IMPORT_NAMES.stream()
                .map(name -> ClassHelper.makeWithoutCaching(name))
                .toList();
        }

        ConfigCompilationUnit(CompilerConfiguration configuration, GroovyClassLoader loader) {
            super(configuration, null, loader);
            super.addPhaseOperation(source -> analyze(source), Phases.CONVERSION);
        }

        private void analyze(SourceUnit source) {
            // initialize script class
            var cn = source.getAST().getClasses().get(0);

            // perform strict syntax checking
            if( source.getSource() instanceof StringReaderSourceWithURI ) {
                var includeResolver = new ResolveIncludeVisitor(source);
                includeResolver.visit();
                for( var error : includeResolver.getErrors() )
                    source.getErrorCollector().addErrorAndContinue(error);
            }
            new ConfigResolveVisitor(source, this, DEFAULT_IMPORTS).visit();
            if( source.getErrorCollector().hasErrors() )
                return;

            // convert to Groovy
            new ConfigToGroovyVisitor(source).visit();
            new PathCompareVisitor(source).visitClass(cn);
            if( stripSecrets )
                new StripSecretsVisitor(source).visitClass(cn);
            if( renderClosureAsString )
                new ClosureToStringVisitor(source).visitClass(cn);
        }

        SourceUnit createSourceUnit(String source, Path path) {
            var readerSource = path != null
                ? new StringReaderSourceWithURI(source, path.toUri(), getConfiguration())
                : new StringReaderSource(source, getConfiguration());
            return new SourceUnit(
                    uniqueClassName(source),
                    readerSource,
                    getConfiguration(),
                    getClassLoader(),
                    getErrorCollector());
        }

        /**
         * Creates a unique name for the config class in order to avoid collision
         * with config DSL
         *
         * @param text
         */
        private String uniqueClassName(String text) {
            var hash = Hashing
                    .sipHash24()
                    .newHasher()
                    .putUnencodedChars(text)
                    .hash();
            return "_nf_config_" + hash;
        }
    }

}
