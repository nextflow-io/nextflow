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

import java.time.Instant

import com.beust.jcommander.IParameterValidator
import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import com.beust.jcommander.ParameterException
import groovy.io.FileType
import groovy.json.JsonOutput
import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.util.logging.Slf4j
import nextflow.config.control.ConfigParser
import nextflow.exception.AbortOperationException
import nextflow.script.control.Compiler
import nextflow.script.control.ScriptParser
import nextflow.util.PathUtils
import org.codehaus.groovy.control.SourceUnit
import org.codehaus.groovy.control.messages.SyntaxErrorMessage
import org.codehaus.groovy.syntax.SyntaxException
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

    @Parameter(
        names = ['-exclude'],
        description = 'File pattern to exclude from linting (can be specified multiple times)'
    )
    List<String> excludePatterns = ['.git', '.nf-test', 'work']

    @Parameter(
        names = ['-output-format'],
        description = 'Output format for lint results [env: NXF_LINT_FORMAT]. Options: full, extended, concise, json',
        validateWith = OutputFormatValidator
    )
    String outputFormat = System.getenv('NXF_LINT_FORMAT') ?: 'full'

    static class OutputFormatValidator implements IParameterValidator {

        private static final List<String> FORMATS = List.of('full', 'extended', 'concise', 'json')

        @Override
        void validate(String name, String value) {
            if( !FORMATS.contains(value) )
                throw new ParameterException("Output format must be one of $FORMATS (found: $value)")
        }
    }

    private ScriptParser scriptParser

    private ConfigParser configParser

    private ErrorListener errorListener

    private LintSummary summary = new LintSummary()

    @Override
    String getName() { 'lint' }

    @Override
    void run() {
        if( !args )
            throw new AbortOperationException("Error: No input files specified")

        scriptParser = new ScriptParser()
        configParser = new ConfigParser()
        errorListener = outputFormat == 'json'
            ? new JsonErrorListener()
            : new AnsiErrorListener(outputFormat)

        errorListener.beforeAll()

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

        errorListener.afterAll(summary)

        // If there were errors, throw an exception to return a non-zero exit code
        if( summary.errors > 0 )
            throw new AbortOperationException()
    }

    void parse(File file) {
        if( PathUtils.isPathExcluded(file.toPath(), excludePatterns) ) {
            log.debug "Skipping excluded file ${file}"
            return
        }

        log.debug "Linting file ${file}"
        errorListener.beforeFile(file)

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
                if( source.getErrorCollector().hasErrors() ) {
                    printErrors(source)
                    summary.filesWithErrors += 1
                }
                else {
                    summary.filesWithoutErrors += 1
                }
            })
    }

    private void printErrors(SourceUnit source) {
        errorListener.beforeErrors()

        final errorMessages = source.getErrorCollector().getErrors()
        for( final message : errorMessages ) {
            if( message instanceof SyntaxErrorMessage ) {
                final cause = message.getCause()
                errorListener.onError(cause, source)
                summary.errors += 1
            }
        }

        errorListener.afterErrors()
    }

}


class LintSummary {
    int errors = 0
    int filesWithErrors = 0
    int filesWithoutErrors = 0
}


interface ErrorListener {
    void beforeAll()
    void beforeFile(File file)
    void beforeErrors()
    void onError(SyntaxException error, SourceUnit source)
    void afterErrors()
    void afterAll(LintSummary summary)
}


@CompileStatic
class AnsiErrorListener implements ErrorListener {
    private String format

    AnsiErrorListener(String format) {
        this.format = format
    }

    @Override
    void beforeAll() {
        final line = ansi().a("Checking Nextflow code..").newline()
        AnsiConsole.out.print(line)
        AnsiConsole.out.flush()
    }

    @Override
    void beforeFile(File file) {
        final line = ansi()
            .cursorUp(1).eraseLine()
            .a(Attribute.INTENSITY_FAINT).a("Checking: ${file.getPath().replaceFirst(/^\.\//, '')}")
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
    void onError(SyntaxException error, SourceUnit source) {
        term.bold().a("${source.getName().replaceFirst(/^\.\//, '')}").reset()
        term.a(":${error.getStartLine()}:${error.getStartColumn()}: ")
        term = highlightString(error.getOriginalMessage(), term)
        if( format != 'concise' ) {
            term.newline()
            term = printCodeBlock(source, error, term)
        }
        term.newline()
    }

    private Ansi highlightString(String str, Ansi term) {
        final matcher = str =~ /^(.*)([`'][^`']+[`'])(.*)$/
        if( matcher.find() ) {
            term.a(matcher.group(1))
                .fg(Color.CYAN).a(matcher.group(2)).fg(Color.DEFAULT)
                .a(matcher.group(3))
        }
        else {
            term.a(str)
        }
        return term
    }

    private Ansi printCodeBlock(SourceUnit source, SyntaxException error, Ansi term) {
        final startLine = error.getStartLine()
        final startColumn = error.getStartColumn()
        final endLine = error.getEndLine()
        final endColumn = error.getEndColumn()
        final lines = getSourceText(source)

        // get context window (up to 5 lines)
        int padding = format == 'extended' ? 2 : 0
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
                // Print line with error in red
                term.a(Attribute.INTENSITY_FAINT).a(line.substring(0, adjStart)).reset()
                term.fg(Ansi.Color.RED).a(line.substring(adjStart, adjEnd)).reset()
                term.a(Attribute.INTENSITY_FAINT).a(line.substring(adjEnd)).reset().newline()

                // Print carets underneath the error range
                String marker = ' ' * adjStart
                String carets = '^' * Math.max(1, adjEnd - adjStart)
                term.a("    | ")
                    .fg(Ansi.Color.RED).bold().a(marker + carets).reset().newline()
            }
            else {
                term.a(Attribute.INTENSITY_FAINT).a(line).reset().newline()
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
        term.fg(Color.DEFAULT).newline()
        AnsiConsole.out.print(term)
        AnsiConsole.out.flush()
    }

    @Override
    void afterAll(LintSummary summary) {
        final emojis = [
            "🔍 📋",
            "🕵🏻‍♂️ 🔎",
            "🔬 👨🏻‍💻",
            "📑 ☑️",
            "🧾 ✔️"
        ]
        final rnd = new Random()
        final term = ansi()
        term.cursorUp(1).eraseLine().cursorUp(1).eraseLine()
        // print extra newline if no code is being shown
        if( format == 'concise' )
            term.newline()
        term.bold().a("Nextflow code checks complete! ${emojis[rnd.nextInt(emojis.size())]}").reset().newline()
        if( summary.filesWithErrors > 0 )
            term.fg(Color.RED).a(" ❌ ${summary.filesWithErrors} file${summary.filesWithErrors==1 ? '' : 's'} had ${summary.errors} error${summary.errors==1 ? '' : 's'}").newline()
        if( summary.filesWithoutErrors > 0 )
            term.fg(Color.BLUE).a(" ✅ ${summary.filesWithoutErrors} file${summary.filesWithoutErrors==1 ? '' : 's'} had no errors").newline()
        if( summary.filesWithErrors == 0 && summary.filesWithoutErrors == 0 )
            term.a(" No files found to process").newline()
        AnsiConsole.out.print(term)
        AnsiConsole.out.flush()
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
    void onError(SyntaxException error, SourceUnit source) {
        errors.add([
            filename: source.getName(),
            startLine: error.getStartLine(),
            startColumn: error.getStartColumn(),
            message: error.getOriginalMessage()
        ])
    }

    @Override
    void afterErrors() {
    }

    @Override
    void afterAll(LintSummary summary) {
        final result = [
            date: Instant.now().toString(),
            summary: summary,
            errors: errors
        ]
        println JsonOutput.prettyPrint(JsonOutput.toJson(result))
    }
}
