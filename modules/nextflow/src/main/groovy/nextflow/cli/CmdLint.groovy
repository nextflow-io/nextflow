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

package nextflow.cli

import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import groovy.io.FileType
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.config.control.ConfigParser
import nextflow.exception.AbortOperationException
import nextflow.script.control.Compiler
import nextflow.script.control.ScriptParser
import org.codehaus.groovy.control.SourceUnit
import org.codehaus.groovy.control.messages.SyntaxErrorMessage
/**
 * CLI sub-command LINT
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
@Parameters(commandDescription = "Check Nextflow scripts and config files for errors")
class CmdLint extends CmdBase {

    @Parameter(description = 'List of paths to lint')
    List<String> args = []

    private ScriptParser scriptParser

    private ConfigParser configParser

    @Override
    String getName() { 'lint' }

    @Override
    void run() {
        if( !args )
            throw new AbortOperationException("No input files specified")

        scriptParser = new ScriptParser()
        configParser = new ConfigParser()

        for( final arg : args ) {
            final file = new File(arg)
            if( file.isFile() ) {
                parse(file)
                continue
            }

            file.eachFileRecurse(FileType.FILES, this.&parse)
        }

        scriptParser.analyze()
        configParser.analyze()
        checkErrors(scriptParser.compiler())
        checkErrors(configParser.compiler())
    }

    void parse(File file) {
        final name = file.getName()
        if( name.endsWith('.nf') )
            scriptParser.parse(file)
        else if( name.endsWith('.config') )
            configParser.parse(file)
    }

    private void checkErrors(Compiler compiler) {
        compiler.getSources()
            .values()
            .stream()
            .sorted(Comparator.comparing((SourceUnit source) -> source.getSource().getURI()))
            .forEach((source) -> {
                if( source.getErrorCollector().hasErrors() )
                    printErrors(source)
            })
    }

    private void printErrors(SourceUnit source) {
        final errorMessages = source.getErrorCollector().getErrors()
        System.err.println()
        for( final message : errorMessages ) {
            if( message instanceof SyntaxErrorMessage ) {
                final cause = message.getCause()
                System.err.println "${source.getName()} at line ${cause.getStartLine()}, column ${cause.getStartColumn()}: ${cause.getOriginalMessage()}"
            }
        }
    }
}
