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
package nextflow.config.control;

import java.io.File;

import groovy.lang.GroovyClassLoader;
import nextflow.config.parser.ConfigParserPluginFactory;
import nextflow.script.control.Compiler;
import nextflow.script.control.LazyErrorCollector;
import nextflow.script.types.Types;
import org.codehaus.groovy.control.CompilerConfiguration;
import org.codehaus.groovy.control.SourceUnit;
import org.codehaus.groovy.control.messages.WarningMessage;

/**
 * Parse and analyze config files without compiling to Groovy.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public class ConfigParser {

    private Compiler compiler;

    public ConfigParser() {
        var config = getConfig();
        var classLoader = new GroovyClassLoader();
        compiler = new Compiler(config, classLoader);
    }

    public Compiler compiler() {
        return compiler;
    }

    public SourceUnit parse(File file) {
        var source = compiler.createSourceUnit(file);
        compiler.addSource(source);
        compiler.compile(source);
        return source;
    }

    public SourceUnit parse(String name, String contents) {
        var source = compiler.createSourceUnit(name, contents);
        compiler.addSource(source);
        compiler.compile(source);
        return source;
    }

    public void analyze() {
        for( var source : compiler.getSources().values() ) {
            var includeResolver = new ResolveIncludeVisitor(source);
            includeResolver.visit();
            for( var error : includeResolver.getErrors() )
                source.getErrorCollector().addErrorAndContinue(error);
            new ConfigResolveVisitor(source, compiler.compilationUnit(), Types.DEFAULT_CONFIG_IMPORTS).visit();
        }
    }

    private static CompilerConfiguration getConfig() {
        var config = new CompilerConfiguration();
        config.setPluginFactory(new ConfigParserPluginFactory());
        config.setWarningLevel(WarningMessage.POSSIBLE_ERRORS);
        return config;
    }

}
