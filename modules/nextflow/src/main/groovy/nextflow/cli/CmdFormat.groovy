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
import java.nio.file.*
import nextflow.config.control.ConfigParser
import nextflow.config.formatter.ConfigFormattingVisitor
import nextflow.exception.AbortOperationException
import nextflow.script.control.ScriptParser
import nextflow.script.formatter.FormattingOptions
import nextflow.script.formatter.ScriptFormattingVisitor
import org.codehaus.groovy.control.SourceUnit
import org.codehaus.groovy.control.messages.SyntaxErrorMessage

import org.fusesource.jansi.Ansi
import org.fusesource.jansi.AnsiConsole
import static nextflow.util.LoggerHelper.isHashLogPrefix
import static org.fusesource.jansi.Ansi.Attribute
import static org.fusesource.jansi.Ansi.Color
import static org.fusesource.jansi.Ansi.ansi
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

    @Parameter(names=['-spaces'], description = 'Number of spaces to indent')
    int spaces

    @Parameter(names = ['-tabs'], description = 'Indent with tabs')
    Boolean tabs

    @Parameter(names = ['-harshil-alignment'], description = 'Use Harshil alignment')
    Boolean harhsilAlignment

    @Parameter(
        names = ['-exclude'],
        description = 'File pattern to exclude from formatter (can be specified multiple times)'
    )
    List<String> excludePatterns = []

    private ScriptParser scriptParser

    private ConfigParser configParser

    private FormattingOptions formattingOptions

    private int numFilesChanged = 0
    private int numFilesUnchanged = 0

    @Override
    String getName() { 'format' }

    @Override
    void run() {
        // Print a newline, as first format update will chomp it
        println()
        if( !args )
            throw new AbortOperationException("Error: No input files specified")

        if( spaces && tabs )
            throw new AbortOperationException("Error: Cannot specify both `-spaces` and `-tabs`")

        if( !spaces && !tabs )
            spaces = 4

        scriptParser = new ScriptParser()
        configParser = new ConfigParser()
        formattingOptions = new FormattingOptions(spaces, !tabs, harhsilAlignment, false)

        for( final arg : args ) {
            final file = new File(arg)
            if( file.isFile() ) {
                format(file)
                continue
            }

            file.eachFileRecurse(FileType.FILES, this.&format)
        }
        final emojis = [
            "🪣 🫧",
            "🧽 ✨",
            "✍️ 💫",
            "🪄 📄",
            "🧼 🎉",
            "🧹 💨"
        ]
        def rnd = new Random()
        final term = ansi().cursorUp(1).eraseLine()
        term.bold().a("Nextflow code formatting complete! ${emojis[rnd.nextInt(emojis.size())]}").reset().newline()
        if(numFilesChanged > 0)
            term.fg(Color.GREEN).a(" ${numFilesChanged} file${numFilesChanged==1 ? '':'s'} reformatted").newline()
        if(numFilesUnchanged > 0)
            term.fg(Color.BLUE).a(" ${numFilesUnchanged} file${numFilesUnchanged==1 ? '':'s'} left unchanged").newline()
        if(numFilesChanged == 0 && numFilesUnchanged == 0)
            term.a(" No files found to process").newline()
        AnsiConsole.out.print(term)
        AnsiConsole.out.flush()
    }

    private boolean shouldExcludeFile(File file) {
        if (!excludePatterns)
            return false

        def path = file.toPath()

        return excludePatterns.any { pattern ->
            if (pattern.contains("*") || pattern.contains("?") || pattern.contains("[")) {
                def matcher = FileSystems.default.getPathMatcher("glob:${pattern}")
                return matcher.matches(path)
            } else {
                def prefix = Paths.get(pattern)
                return path.getNameCount() >= prefix.getNameCount() &&
                    path.subpath(0, prefix.getNameCount()) == prefix
            }
        }
    }

    void format(File file) {
        final name = file.getName()
        // Skip excluded files
        if (shouldExcludeFile(file)) {
            log.debug "Skipping excluded file ${file}"
            return
        }

        if( name.endsWith('.nf') )
            formatScript(file)
        else if( name.endsWith('.config') )
            formatConfig(file)
    }

    private void formatScript(File file) {
        final source = scriptParser.parse(file)
        if( source.getErrorCollector().hasErrors() ) {
            printErrors(source)
        }
        else {
            log.debug "Formatting script ${file}"
            printStatus(file)
            final formatter = new ScriptFormattingVisitor(source, formattingOptions)
            formatter.visit()
            // Check if anything changed)
            if(file.text != formatter.toString()){
                numFilesChanged += 1
                file.text = formatter.toString()
            } else {
                numFilesUnchanged += 1
            }
        }
    }

    private void formatConfig(File file) {
        final source = configParser.parse(file)
        if( source.getErrorCollector().hasErrors() ) {
            printErrors(source)
        }
        else {
            log.debug "Formatting config ${file.getPath().replaceFirst(/^\.\//, '')}"
            printStatus(file)
            final formatter = new ConfigFormattingVisitor(source, formattingOptions)
            formatter.visit()
            // Check if anything changed)
            if(file.text != formatter.toString()){
                numFilesChanged += 1
                file.text = formatter.toString()
            } else {
                numFilesUnchanged += 1
            }
        }
    }

    private void printStatus(File file) {
        final str = ansi().cursorUp(1).eraseLine()
        str.a(Attribute.INTENSITY_FAINT).a("Formatting: ${file.getPath().replaceFirst(/^\.\//, '')}")
        str.reset().newline().toString()
        AnsiConsole.out.print(str)
        AnsiConsole.out.flush()
    }

    private Ansi highlightString(String str, Ansi term) {
        def matcher = str =~ /^(.*)([`'][^`']+[`'])(.*)$/
        if (matcher.find()) {
            term.a(matcher.group(1))
                .fg(Color.CYAN).a(matcher.group(2)).fg(Color.DEFAULT)
                .a(matcher.group(3))
        } else {
            term.a(str)
        }
        return term
    }

    private void printErrors(SourceUnit source) {
        final errorMessages = source.getErrorCollector().getErrors()
        def term = ansi().cursorUp(1).eraseLine()
        for( final message : errorMessages ) {
            if( message instanceof SyntaxErrorMessage ) {
                final cause = message.getCause()
                term.fg(Color.RED).bold().a("error").fg(Color.DEFAULT).a(": ")
                term.a("Failed to parse ${source.getName().replaceFirst(/^\.\//, '')}").a(Attribute.INTENSITY_BOLD_OFF)
                term.a(":${cause.getStartLine()}:${cause.getStartColumn()}: ")
                term = highlightString(cause.getOriginalMessage(), term)
                term.newline()
            }
        }
        // Double newline as next status update will chomp back one
        term.fg(Color.DEFAULT).newline()
        AnsiConsole.out.print(term.toString())
        AnsiConsole.out.flush()
    }
}
