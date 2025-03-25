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
import org.fusesource.jansi.Ansi
import org.fusesource.jansi.AnsiConsole
import static org.fusesource.jansi.Ansi.Attribute
import static org.fusesource.jansi.Ansi.Color
import static org.fusesource.jansi.Ansi.ansi
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

    private int numErrors = 0
    private int numFilesWithErrors = 0
    private int numFilesWithoutErrors = 0

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
        printStatus(file)
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
                if( source.getErrorCollector().hasErrors() ){
                    printErrors(source)
                    numFilesWithErrors += 1
                } else {
                    numFilesWithoutErrors += 1
                }
            })
    }

    private void printStatus(File file) {
        final str = ansi().cursorUp(1).eraseLine().a(Attribute.INTENSITY_FAINT).a("Formatting: ${file}").reset().newline().toString()
        AnsiConsole.out.print(str)
        AnsiConsole.out.flush()
    }

    private void printErrors(SourceUnit source) {
        final errorMessages = source.getErrorCollector().getErrors()
        final term = ansi().cursorUp(1).eraseLine()
        for( final message : errorMessages ) {
            if( message instanceof SyntaxErrorMessage ) {
                numErrors += 1
                final cause = message.getCause()
                term.fg(Color.RED).a(Attribute.INTENSITY_BOLD).a("error").fg(Color.DEFAULT).a(": ")
                term.a("${source.getName()}").a(Attribute.INTENSITY_BOLD_OFF)
                term.a(":${cause.getStartLine()}:${cause.getStartColumn()}: ${cause.getOriginalMessage()}")
                term.newline()
            }
        }
        // Double newline as next status update will chomp back one
        term.fg(Color.DEFAULT).newline()
        AnsiConsole.out.print(term)
        AnsiConsole.out.flush()
    }
}
