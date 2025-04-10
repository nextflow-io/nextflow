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

import java.nio.file.Path

import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import groovy.io.FileType
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.config.control.ConfigParser
import nextflow.config.formatter.ConfigFormattingVisitor
import nextflow.exception.AbortOperationException
import nextflow.script.control.ScriptParser
import nextflow.script.formatter.FormattingOptions
import nextflow.script.formatter.ScriptFormattingVisitor
import nextflow.util.PathUtils
import org.codehaus.groovy.control.SourceUnit
import org.codehaus.groovy.control.messages.SyntaxErrorMessage
import org.fusesource.jansi.Ansi
import org.fusesource.jansi.AnsiConsole
/**
 * CLI sub-command FORMAT
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
@Parameters(commandDescription = "Format Nextflow scripts and config files")
class CmdFormat extends CmdBase {

    @Parameter(description = 'List of paths to format')
    List<String> args = []

    @Parameter(
        names = ['-exclude'],
        description = 'File pattern to exclude from formatter (can be specified multiple times)'
    )
    List<String> excludePatterns = ['.git', '.nf-test', 'work']

    @Parameter(names = ['-harshil-alignment'], description = 'Use Harshil alignment')
    Boolean harhsilAlignment

    @Parameter(names = ['-sort-declarations'], description = 'Sort script declarations in Nextflow scripts')
    Boolean sortDeclarations

    @Parameter(names=['-spaces'], description = 'Number of spaces to indent')
    int spaces

    @Parameter(names = ['-tabs'], description = 'Indent with tabs')
    Boolean tabs

    private ScriptParser scriptParser

    private ConfigParser configParser

    private FormattingOptions formattingOptions

    private int filesChanged = 0
    private int filesUnchanged = 0

    @Override
    String getName() { 'format' }

    @Override
    void run() {
        if( !args )
            throw new AbortOperationException("Error: No input files were specified")

        if( spaces && tabs )
            throw new AbortOperationException("Error: Cannot specify both `-spaces` and `-tabs`")

        if( !spaces && !tabs )
            spaces = 4

        scriptParser = new ScriptParser()
        configParser = new ConfigParser()
        formattingOptions = new FormattingOptions(spaces, !tabs, harhsilAlignment, false, sortDeclarations)

        // print extra newline since first file status will chomp it
        println()

        for( final arg : args ) {
            PathUtils.visitFiles(
                Path.of(arg),
                (path) -> !PathUtils.isExcluded(path, excludePatterns),
                (path) -> format(path.toFile()))
        }

        final term = ansi().cursorUp(1).eraseLine()
        term.bold().a("Nextflow code formatting complete!").reset().newline()
        if( filesChanged > 0 )
            term.fg(Ansi.Color.GREEN).a(" ${filesChanged} file${filesChanged==1 ? '':'s'} reformatted").newline()
        if( filesUnchanged > 0 )
            term.fg(Ansi.Color.BLUE).a(" ${filesUnchanged} file${filesUnchanged==1 ? '':'s'} left unchanged").newline()
        if( filesChanged == 0 && filesUnchanged == 0 )
            term.a(" No files found to process").newline()
        AnsiConsole.out.print(term)
        AnsiConsole.out.flush()
    }

    private Ansi ansi() {
        final ansi = Ansi.ansi()
        ansi.setEnabled(launcher.options.ansiLog)
        return ansi
    }

    private void format(File file) {
        final name = file.getName()
        if( name.endsWith('.nf') )
            formatScript(file)
        else if( name.endsWith('.config') )
            formatConfig(file)
    }

    private void formatScript(File file) {
        final source = scriptParser.parse(file)
        if( source.getErrorCollector().hasErrors() ) {
            printErrors(source)
            return
        }

        log.debug "Formatting script ${file}"
        printStatus(file)

        final formatter = new ScriptFormattingVisitor(source, formattingOptions)
        formatter.visit()
        final result = formatter.toString()
        if( file.text != result ) {
            filesChanged += 1
            file.text = result
        }
        else {
            filesUnchanged += 1
        }
    }

    private void formatConfig(File file) {
        final source = configParser.parse(file)
        if( source.getErrorCollector().hasErrors() ) {
            printErrors(source)
            return
        }

        log.debug "Formatting config ${file}"
        printStatus(file)

        final formatter = new ConfigFormattingVisitor(source, formattingOptions)
        formatter.visit()
        final result = formatter.toString()
        if( file.text != result ) {
            filesChanged += 1
            file.text = result
        }
        else {
            filesUnchanged += 1
        }
    }

    private void printStatus(File file) {
        final line = ansi()
            .cursorUp(1).eraseLine()
            .a(Ansi.Attribute.INTENSITY_FAINT).a("Formatting: ${file}")
            .reset().newline().toString()
        AnsiConsole.out.print(line)
        AnsiConsole.out.flush()
    }

    private void printErrors(SourceUnit source) {
        Ansi term = ansi().cursorUp(1).eraseLine()

        final errorMessages = source.getErrorCollector().getErrors()
        for( final message : errorMessages ) {
            if( message instanceof SyntaxErrorMessage ) {
                final cause = message.getCause()
                term.fg(Ansi.Color.RED).bold().a("error").fg(Ansi.Color.DEFAULT).a(": ")
                term.a("Failed to parse ${source.getName()}").a(Ansi.Attribute.INTENSITY_BOLD_OFF)
                term.a(":${cause.getStartLine()}:${cause.getStartColumn()}: ")
                term = highlightString(cause.getOriginalMessage(), term)
                term.newline()
            }
        }

        // print extra newline since next file status will chomp back one
        term.fg(Ansi.Color.DEFAULT).newline()
        AnsiConsole.out.print(term.toString())
        AnsiConsole.out.flush()
    }

    private Ansi highlightString(String str, Ansi term) {
        final matcher = str =~ /^(.*)([`'][^`']+[`'])(.*)$/
        if( matcher.find() ) {
            term.a(matcher.group(1))
                .fg(Ansi.Color.CYAN).a(matcher.group(2)).fg(Ansi.Color.DEFAULT)
                .a(matcher.group(3))
        }
        else {
            term.a(str)
        }
        return term
    }
}
