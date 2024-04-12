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
package nextflow.script.v2

import nextflow.exception.ScriptCompilationException
import org.codehaus.groovy.ast.ModuleNode
import org.codehaus.groovy.control.ParserPlugin
import org.codehaus.groovy.control.SourceUnit
import org.codehaus.groovy.control.io.StringReaderSource
import org.codehaus.groovy.runtime.IOGroovyMethods
import org.codehaus.groovy.syntax.Reduction

/**
 * Parser plugin for the Nextflow parser.
 */
class ScriptParserPlugin implements ParserPlugin {

    @Override
    Reduction parseCST(SourceUnit sourceUnit, Reader reader) {
        if (!sourceUnit.getSource().canReopenSource()) {
            try {
                sourceUnit.setSource(new StringReaderSource(
                        IOGroovyMethods.getText(reader),
                        sourceUnit.getConfiguration()
                ))
            } catch (IOException e) {
                throw new ScriptCompilationException("Failed to create StringReaderSource", e)
            }
        }
        return null
    }

    @Override
    ModuleNode buildAST(SourceUnit sourceUnit, ClassLoader classLoader, Reduction cst) {
        return new AstBuilder(sourceUnit).buildAST()
    }
}
