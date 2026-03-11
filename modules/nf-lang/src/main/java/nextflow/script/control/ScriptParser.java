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
import java.util.ArrayList;
import java.util.Collections;

import groovy.lang.GroovyClassLoader;
import nextflow.script.parser.ScriptParserPluginFactory;
import nextflow.script.types.Types;
import org.codehaus.groovy.control.CompilerConfiguration;
import org.codehaus.groovy.control.SourceUnit;
import org.codehaus.groovy.control.messages.WarningMessage;

/**
 * Parse and analyze scripts without compiling to Groovy.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public class ScriptParser {

    private Compiler compiler;

    public ScriptParser() {
        var config = getConfig();
        var classLoader = new GroovyClassLoader();
        compiler = new Compiler(config, classLoader);
    }

    public Compiler compiler() {
        return compiler;
    }

    public SourceUnit parse(File file) {
        var source = compiler.createSourceUnit(file);
        parse0(source);
        return source;
    }

    public SourceUnit parse(String name, String contents) {
        var source = compiler.createSourceUnit(name, contents);
        parse0(source);
        return source;
    }

    private void parse0(SourceUnit source) {
        var uri = source.getSource().getURI();
        if( compiler.getSource(uri) != null )
            return;
        compiler.addSource(source);
        compiler.compile(source);
    }

    public void analyze() {
        var sources = new ArrayList<>(compiler.getSources().values());
        for( var source : sources ) {
            new ModuleResolver(compiler()).resolve(source, (uri) -> compiler.createSourceUnit(new File(uri)));
        }

        for( var source : compiler.getSources().values() ) {
            var includeResolver = new ResolveIncludeVisitor(source, compiler);
            includeResolver.visit();
            for( var error : includeResolver.getErrors() )
                source.getErrorCollector().addErrorAndContinue(error);
            new ScriptResolveVisitor(source, compiler.compilationUnit(), Types.DEFAULT_SCRIPT_IMPORTS, Collections.emptyList()).visit();
            if( source.getErrorCollector().hasErrors() )
                continue;
            new TypeCheckingVisitor(source).visit();
        }
    }

    private static CompilerConfiguration getConfig() {
        var config = new CompilerConfiguration();
        config.setPluginFactory(new ScriptParserPluginFactory());
        config.setWarningLevel(WarningMessage.POSSIBLE_ERRORS);
        return config;
    }

}
