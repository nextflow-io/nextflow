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
package nextflow.config.parser;

import org.codehaus.groovy.GroovyBugError;
import org.codehaus.groovy.ast.ModuleNode;
import org.codehaus.groovy.control.ParserPlugin;
import org.codehaus.groovy.control.SourceUnit;
import org.codehaus.groovy.control.io.StringReaderSource;
import org.codehaus.groovy.runtime.IOGroovyMethods;
import org.codehaus.groovy.syntax.Reduction;

import java.io.IOException;
import java.io.Reader;

/**
 * Parser plugin for the Nextflow config parser.
 */
public class ConfigParserPlugin implements ParserPlugin {

    @Override
    public Reduction parseCST(SourceUnit sourceUnit, Reader reader) {
        if (!sourceUnit.getSource().canReopenSource()) {
            try {
                sourceUnit.setSource(new StringReaderSource(
                        IOGroovyMethods.getText(reader),
                        sourceUnit.getConfiguration()
                ));
            } catch (IOException e) {
                throw new GroovyBugError("Failed to create StringReaderSource", e);
            }
        }
        return null;
    }

    @Override
    public ModuleNode buildAST(SourceUnit sourceUnit, ClassLoader classLoader, Reduction cst) {
        return new ConfigAstBuilder(sourceUnit).buildAST();
    }
}
