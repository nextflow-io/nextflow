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
import java.time.Instant

import com.beust.jcommander.IParameterValidator
import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import com.beust.jcommander.ParameterException
import groovy.json.JsonOutput
import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.util.logging.Slf4j
import nextflow.config.control.ConfigParser
import nextflow.config.formatter.ConfigFormattingVisitor
import nextflow.exception.AbortOperationException
import nextflow.script.control.Compiler
import nextflow.script.control.ParanoidWarning
import nextflow.script.control.ScriptParser
import nextflow.script.formatter.FormattingOptions
import nextflow.script.formatter.ScriptFormattingVisitor
import nextflow.util.PathUtils
import org.codehaus.groovy.control.SourceUnit
import org.codehaus.groovy.control.messages.SyntaxErrorMessage
import org.codehaus.groovy.control.messages.WarningMessage
import org.codehaus.groovy.syntax.SyntaxException
import org.fusesource.jansi.Ansi
import org.fusesource.jansi.AnsiConsole
/**
 * CLI sub-command LINT
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
@Parameters(commandDescription = "Lint Nextflow scripts and config files")
class CmdLint extends CmdBase {

    @Parameter(description = 'List of paths to lint')
    List<String> args = []

    @Parameter(
        names = ['-exclude'],
        description = 'File pattern to exclude from error checking (can be specified multiple times)'
    )
    List<String> excludePatterns = ['.git', '.nf-test', 'work']

    @Parameter(
        names = ['-o', '-output'],
        description = 'Output mode for reporting errors: full, extended, concise, json',
        validateWith = OutputModeValidator
    )
    String outputMode = 'full'

    static class OutputModeValidator implements IParameterValidator {

        private static final List<String> MODES = List.of('full', 'extended', 'concise', 'json')

        @Override
        void validate(String name, String value) {
            if( !MODES.contains(value) )
                throw new ParameterException("Output mode must be one of $MODES (found: $value)")
        }
    }

    @Parameter(names = ['-format'], description = 'Format scripts and config files that have no errors')
    boolean formatting

    @Parameter(names = ['-harshil-alignment'], description = 'Use Harshil alignment')
    boolean harhsilAlignment

    @Parameter(names = ['-sort-declarations'], description = 'Sort script declarations in Nextflow scripts')
    boolean sortDeclarations

    @Parameter(names=['-spaces'], description = 'Number of spaces to indent')
    int spaces

    @Parameter(names = ['-tabs'], description = 'Indent with tabs')
    boolean tabs

    private ScriptParser scriptParser

    private ConfigParser configParser

    private ErrorListener errorListener

    private FormattingOptions formattingOptions

    private ErrorSummary summary = new ErrorSummary()

    @Override
    String getName() { 'lint' }

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
        errorListener = outputMode == 'json'
            ? new JsonErrorListener()
            : new StandardErrorListener(outputMode, launcher.options.ansiLog)
        formattingOptions = new FormattingOptions(spaces, !tabs, harhsilAlignment, false, sortDeclarations)

        errorListener.beforeAll()

        // collect files to lint
        final List<File> files = []

        for( final arg : args ) {
            PathUtils.visitFiles(
                Path.of(arg),
                (path) -> !PathUtils.isExcluded(path, excludePatterns),
                (path) -> files.add(path.toFile()))
        }

        // parse and analyze files
        for( final file : files )
            parse(file)

        scriptParser.analyze()
        configParser.analyze()

        // report errors
        checkErrors(scriptParser.compiler())
        checkErrors(configParser.compiler())

        // format files if specified
        if( formatting ) {
            for( final file : files )
                format(file)
        }

        // print summary
        errorListener.afterAll(summary)

        // If there were errors, throw an exception to return a non-zero exit code
        if( summary.errors > 0 )
            throw new AbortOperationException()
    }

    private void parse(File file) {
        final name = file.getName()
        if( name.endsWith('.nf') )
            parseScript(file)
        else if( name.endsWith('.config') )
            parseConfig(file)
    }

    private void parseScript(File file) {
        log.debug "Linting script ${file}"
        errorListener.beforeFile(file)
        scriptParser.parse(file)
    }

    private void parseConfig(File file) {
        log.debug "Linting config ${file}"
        errorListener.beforeFile(file)
        configParser.parse(file)
    }

    private void checkErrors(Compiler compiler) {
        compiler.getSources()
            .values()
            .stream()
            .sorted(Comparator.comparing((SourceUnit source) -> source.getSource().getURI()))
            .forEach((source) -> {
                final errorCollector = source.getErrorCollector()
                if( errorCollector.hasErrors() || errorCollector.hasWarnings() )
                    printErrors(source)
                if( errorCollector.hasErrors() )
                    summary.filesWithErrors += 1
                else
                    summary.filesWithoutErrors += 1
            })
    }

    private void printErrors(SourceUnit source) {
        errorListener.beforeErrors()

        final errorMessages = source.getErrorCollector().getErrors()
        for( final message : errorMessages ) {
            if( message instanceof SyntaxErrorMessage ) {
                final cause = message.getCause()
                errorListener.onError(cause, source.getName(), source)
                summary.errors += 1
            }
        }

        final warningMessages = source.getErrorCollector().getWarnings()
        for( final warning : warningMessages ) {
            if( warning instanceof ParanoidWarning )
                continue
            errorListener.onWarning(warning, source.getName(), source)
        }

        errorListener.afterErrors()
    }

    private void format(File file) {
        final name = file.getName()
        final result =
            name.endsWith('.nf') ? formatScript(file) :
            name.endsWith('.config') ? formatConfig(file) :
            null

        if( result != null && file.text != result ) {
            summary.filesFormatted += 1
            file.text = result
        }
    }

    private String formatScript(File file) {
        final source = scriptParser.compiler().getSource(file.toURI())
        if( source.getErrorCollector().hasErrors() ) {
            printErrors(source)
            return null
        }

        log.debug "Formatting script ${file}"
        errorListener.beforeFormat(file)

        final formatter = new ScriptFormattingVisitor(source, formattingOptions)
        formatter.visit()
        return formatter.toString()
    }

    private String formatConfig(File file) {
        final source = configParser.compiler().getSource(file.toURI())
        if( source.getErrorCollector().hasErrors() ) {
            printErrors(source)
            return null
        }

        log.debug "Formatting config ${file}"
        errorListener.beforeFormat(file)

        final formatter = new ConfigFormattingVisitor(source, formattingOptions)
        formatter.visit()
        return formatter.toString()
    }

}


class ErrorSummary {
    int errors = 0
    int filesWithErrors = 0
    int filesWithoutErrors = 0
    int filesFormatted = 0
}


interface ErrorListener {
    void beforeAll()
    void beforeFile(File file)
    void beforeErrors()
    void onError(SyntaxException error, String filename, SourceUnit source)
    void onWarning(WarningMessage warning, String filename, SourceUnit source)
    void afterErrors()
    void beforeFormat(File file)
    void afterAll(ErrorSummary summary)
}


@CompileStatic
class StandardErrorListener implements ErrorListener {
    private String mode
    private boolean ansiLog

    StandardErrorListener(String mode, boolean ansiLog) {
        this.mode = mode
        this.ansiLog = ansiLog
    }

    private Ansi ansi() {
        final ansi = Ansi.ansi()
        ansi.setEnabled(ansiLog)
        return ansi
    }

    @Override
    void beforeAll() {
        final line = ansi().a("Linting Nextflow code..").newline()
        AnsiConsole.out.print(line)
        AnsiConsole.out.flush()
    }

    @Override
    void beforeFile(File file) {
        final line = ansi()
            .cursorUp(1).eraseLine()
            .a(Ansi.Attribute.INTENSITY_FAINT).a("Linting: ${file}")
            .reset().newline().toString()
        AnsiConsole.out.print(line)
        AnsiConsole.out.flush()
    }

    private Ansi term

    @Override
    void beforeErrors() {
        term = ansi().cursorUp(1).eraseLine()
    }

    @Override
    void onError(SyntaxException error, String filename, SourceUnit source) {
        term.bold().a(filename).reset()
        term.a(":${error.getStartLine()}:${error.getStartColumn()}: ")
        term = highlightString(error.getOriginalMessage(), term)
        if( mode != 'concise' ) {
            term.newline()
            term = printCodeBlock(source, Range.of(error), term, Ansi.Color.RED)
        }
        term.newline()
    }

    @Override
    void onWarning(WarningMessage warning, String filename, SourceUnit source) {
        final token = warning.getContext().getRoot()
        term.bold().a(filename).reset()
        term.a(":${token.getStartLine()}:${token.getStartColumn()}: ")
        term.fg(Ansi.Color.YELLOW).a(warning.getMessage()).fg(Ansi.Color.DEFAULT)
        if( mode != 'concise' ) {
            term.newline()
            term = printCodeBlock(source, Range.of(warning), term, Ansi.Color.YELLOW)
        }
        term.newline()
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

    private Ansi printCodeBlock(SourceUnit source, Range range, Ansi term, Ansi.Color color) {
        final startLine = range.startLine()
        final startColumn = range.startColumn()
        final endLine = range.endLine()
        final endColumn = range.endColumn()
        final lines = getSourceText(source)

        // get context window (up to 5 lines)
        int padding = mode == 'extended' ? 2 : 0
        int fromLine = Math.max(1, startLine - padding)
        int toLine = Math.min(lines.size(), endLine + padding)
        if( toLine - fromLine + 1 > 5 ) {
            if( startLine <= 3 ) {
                toLine = fromLine + 4
            }
            else if( endLine >= lines.size() - 2 ) {
                fromLine = toLine - 4
            }
            else {
                fromLine = startLine - 2
                toLine = startLine + 2
            }
        }

        for( int i = fromLine; i <= toLine; i++ ) {
            String fullLine = lines[i - 1]
            int start = (i == startLine) ? startColumn - 1 : 0
            int end = (i == endLine) ? endColumn - 1 : fullLine.length()

            // Truncate to max 70 characters
            int maxLen = 70
            int lineLen = fullLine.length()
            int windowStart = 0
            if( lineLen > maxLen ) {
                if( start < maxLen - 10 )
                    windowStart = 0
                else if( end > lineLen - 10 )
                    windowStart = lineLen - maxLen
                else
                    windowStart = start - 30
            }

            String line = fullLine.substring(windowStart, Math.min(lineLen, windowStart + maxLen))
            int adjStart = Math.max(0, start - windowStart)
            int adjEnd = Math.max(adjStart + 1, Math.min(end - windowStart, line.length()))

            // Line number
            term.fg(Ansi.Color.BLUE).a(String.format("%3d | ", i)).reset()

            if( i == startLine ) {
                // Print line with range highlighted
                term.a(Ansi.Attribute.INTENSITY_FAINT).a(line.substring(0, adjStart)).reset()
                term.fg(color).a(line.substring(adjStart, adjEnd)).reset()
                term.a(Ansi.Attribute.INTENSITY_FAINT).a(line.substring(adjEnd)).reset().newline()

                // Print carets underneath the range
                String marker = ' ' * adjStart
                String carets = '^' * Math.max(1, adjEnd - adjStart)
                term.a("    | ")
                    .fg(color).bold().a(marker + carets).reset().newline()
            }
            else {
                term.a(Ansi.Attribute.INTENSITY_FAINT).a(line).reset().newline()
            }
        }

        return term
    }

    @Memoized
    private List<String> getSourceText(SourceUnit source) {
        return source.getSource().getReader().readLines()
    }

    @Override
    void afterErrors() {
        // print extra newline since next file status will chomp back one
        term.fg(Ansi.Color.DEFAULT).newline()
        AnsiConsole.out.print(term)
        AnsiConsole.out.flush()
    }

    @Override
    void beforeFormat(File file) {
        final line = ansi()
            .cursorUp(1).eraseLine()
            .a(Ansi.Attribute.INTENSITY_FAINT).a("Formatting: ${file}")
            .reset().newline().toString()
        AnsiConsole.out.print(line)
        AnsiConsole.out.flush()
    }

    @Override
    void afterAll(ErrorSummary summary) {
        final term = ansi()
        term.cursorUp(1).eraseLine().cursorUp(1).eraseLine()
        // print extra newline if no code is being shown
        if( mode == 'concise' )
            term.newline()
        term.bold().a("Nextflow linting complete!").reset().newline()
        if( summary.filesWithErrors > 0 ) {
            term.fg(Ansi.Color.RED).a(" ❌ ${summary.filesWithErrors} file${summary.filesWithErrors==1 ? '' : 's'} had ${summary.errors} error${summary.errors==1 ? '' : 's'}").newline()
        }
        if( summary.filesWithoutErrors > 0 ) {
            term.fg(Ansi.Color.GREEN).a(" ✅ ${summary.filesWithoutErrors} file${summary.filesWithoutErrors==1 ? '' : 's'} had no errors")
            if( summary.filesFormatted > 0 )
                term.fg(Ansi.Color.BLUE).a(" (${summary.filesFormatted} formatted)")
            term.newline()
        }
        if( summary.filesWithErrors == 0 && summary.filesWithoutErrors == 0 ) {
            term.a(" No files found to process").newline()
        }
        AnsiConsole.out.print(term)
        AnsiConsole.out.flush()
    }

    private static record Range(
        int startLine,
        int startColumn,
        int endLine,
        int endColumn
    ) {
        public static Range of(SyntaxException error) {
            return new Range(
                error.getStartLine(),
                error.getStartColumn(),
                error.getEndLine(),
                error.getEndColumn(),
            )
        }

        public static Range of(WarningMessage warning) {
            final token = warning.getContext().getRoot()
            return new Range(
                token.getStartLine(),
                token.getStartColumn(),
                token.getStartLine(),
                token.getStartColumn() + token.getText().length(),
            )
        }
    }
}


@CompileStatic
class JsonErrorListener implements ErrorListener {

    private List<Map> errors = []

    @Override
    void beforeAll() {
    }

    @Override
    void beforeFile(File file) {
    }

    @Override
    void beforeErrors() {
    }

    @Override
    void onError(SyntaxException error, String filename, SourceUnit source) {
        errors.add([
            filename: filename,
            startLine: error.getStartLine(),
            startColumn: error.getStartColumn(),
            message: error.getOriginalMessage()
        ])
    }

    @Override
    void onWarning(WarningMessage warning, String filename, SourceUnit source) {
    }

    @Override
    void afterErrors() {
    }

    @Override
    void beforeFormat(File file) {
    }

    @Override
    void afterAll(ErrorSummary summary) {
        final result = [
            date: Instant.now().toString(),
            summary: summary,
            errors: errors
        ]
        println JsonOutput.prettyPrint(JsonOutput.toJson(result))
    }
}
