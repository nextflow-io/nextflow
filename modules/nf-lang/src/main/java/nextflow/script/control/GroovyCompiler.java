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

import java.io.File;
import java.nio.file.Path;
import java.util.Collections;
import java.util.List;

import groovy.lang.GroovyClassLoader;
import org.codehaus.groovy.GroovyBugError;
import org.codehaus.groovy.ast.ClassNode;
import org.codehaus.groovy.control.CompilationUnit;
import org.codehaus.groovy.control.CompilerConfiguration;
import org.codehaus.groovy.control.Phases;
import org.codehaus.groovy.control.SourceUnit;

/**
 * Load Groovy classes from the `lib` directory.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public class GroovyCompiler {

    public static List<ClassNode> compile(SourceUnit su) {
        // create compilation unit
        var config = new CompilerConfiguration();
        config.getOptimizationOptions().put(CompilerConfiguration.GROOVYDOC, true);
        var classLoader = new GroovyClassLoader();
        var compilationUnit = new CompilationUnit(config, null, classLoader);

        // create source units (or restore from cache)
        var uri = su.getSource().getURI();
        var sourceUnit = new SourceUnit(
                new File(uri),
                config,
                classLoader,
                new LazyErrorCollector(config));
        compilationUnit.addSource(sourceUnit);

        // compile source files
        try {
            compilationUnit.compile(Phases.CANONICALIZATION);
        }
        catch( GroovyBugError | Exception e ) {
            // ignore
        }

        // collect compiled classes
        var moduleNode = sourceUnit.getAST();
        if( moduleNode == null )
            return Collections.emptyList();

        return moduleNode.getClasses();
    }

}
