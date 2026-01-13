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
import java.time.ZoneOffset
import java.time.format.DateTimeFormatter

import com.beust.jcommander.IParameterValidator
import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import com.beust.jcommander.ParameterException
import groovy.json.JsonOutput
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.BuildInfo
import nextflow.config.control.ConfigParser
import nextflow.config.formatter.ConfigFormattingVisitor
import nextflow.exception.AbortOperationException
import nextflow.script.control.Compiler
import nextflow.script.control.ParanoidWarning
import nextflow.script.control.ScriptParser
import nextflow.script.formatter.FormattingOptions
import nextflow.script.formatter.ScriptFormattingVisitor
import nextflow.script.parser.v2.ErrorListener
import nextflow.script.parser.v2.ErrorSummary
import nextflow.script.parser.v2.StandardErrorListener
import nextflow.util.PathUtils
import org.codehaus.groovy.control.SourceUnit
import org.codehaus.groovy.control.messages.SyntaxErrorMessage
import org.codehaus.groovy.control.messages.WarningMessage
import org.codehaus.groovy.syntax.SyntaxException
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
    List<String> excludePatterns = ['.git', '.lineage', '.nf-test', '.nextflow', 'work']

    @Parameter(
        names = ['-o', '-output'],
        description = 'Output mode for reporting errors: full, extended, concise, json, markdown',
        validateWith = OutputModeValidator
    )
    String outputMode = 'full'

    static class OutputModeValidator implements IParameterValidator {

        private static final List<String> MODES = List.of('full', 'extended', 'concise', 'json', 'markdown')

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
            : outputMode == 'markdown'
            ? new MarkdownErrorListener()
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
                final hasWarnings = (errorCollector.getWarnings() ?: []).stream()
                    .anyMatch(warning -> warning !instanceof ParanoidWarning)
                if( errorCollector.hasErrors() || hasWarnings )
                    printErrors(source)
                if( errorCollector.hasErrors() )
                    summary.filesWithErrors += 1
                else
                    summary.filesWithoutErrors += 1
                if( hasWarnings )
                    summary.filesWithWarnings += 1
                else
                    summary.filesWithoutWarnings += 1
            })
    }

    private void printErrors(SourceUnit source) {
        errorListener.beforeErrors()

        final errors = source.getErrorCollector().getErrors() ?: []
        errors.stream()
            .filter(message -> message instanceof SyntaxErrorMessage)
            .map(message -> ((SyntaxErrorMessage) message).getCause())
            .sorted(ERROR_COMPARATOR)
            .forEach((cause) -> {
                errorListener.onError(cause, source.getName(), source)
                summary.errors += 1
            })

        final warnings = source.getErrorCollector().getWarnings() ?: []
        warnings.stream()
            .filter(warning -> warning !instanceof ParanoidWarning)
            .sorted(WARNING_COMPARATOR)
            .forEach((warning) -> {
                errorListener.onWarning(warning, source.getName(), source)
                summary.warnings += 1
            })

        errorListener.afterErrors()
    }

    private static final Comparator<SyntaxException> ERROR_COMPARATOR = (SyntaxException a, SyntaxException b) -> {
        return a.getStartLine() != b.getStartLine()
            ? a.getStartLine() - b.getStartLine()
            : a.getStartColumn() - b.getStartColumn()
    }

    private static final Comparator<WarningMessage> WARNING_COMPARATOR = (WarningMessage w1, WarningMessage w2) -> {
        final a = w1.getContext()
        final b = w2.getContext()
        return a.getStartLine() != b.getStartLine()
            ? a.getStartLine() - b.getStartLine()
            : a.getStartColumn() - b.getStartColumn()
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


@CompileStatic
class JsonErrorListener implements ErrorListener {

    private List<Map> errors = []

    private List<Map> warnings = []

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
        final token = warning.getContext().getRoot()
        warnings.add([
            filename: filename,
            startLine: token.getStartLine(),
            startColumn: token.getStartColumn(),
            message: warning.getMessage()
        ])
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
            errors: errors,
            warnings: warnings
        ]
        println JsonOutput.prettyPrint(JsonOutput.toJson(result))
    }
}


@CompileStatic
class MarkdownErrorListener implements ErrorListener {

    private static class LintEntry {
        String filename
        int startLine
        int startColumn
        int endLine
        int endColumn
        String message
        SourceUnit source
    }

    private List<LintEntry> errors = []

    private List<LintEntry> warnings = []

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
        errors.add(new LintEntry(
            filename: filename,
            startLine: error.getStartLine(),
            startColumn: error.getStartColumn(),
            endLine: error.getEndLine(),
            endColumn: error.getEndColumn(),
            message: error.getOriginalMessage(),
            source: source
        ))
    }

    @Override
    void onWarning(WarningMessage warning, String filename, SourceUnit source) {
        final token = warning.getContext().getRoot()
        warnings.add(new LintEntry(
            filename: filename,
            startLine: token.getStartLine(),
            startColumn: token.getStartColumn(),
            endLine: token.getStartLine(),
            endColumn: token.getStartColumn() + token.getText().length(),
            message: warning.getMessage(),
            source: source
        ))
    }

    @Override
    void afterErrors() {
    }

    @Override
    void beforeFormat(File file) {
    }

    @Override
    void afterAll(ErrorSummary summary) {
        final sb = new StringBuilder()

        // Header
        sb.append('# Nextflow lint results\n\n')

        // Metadata
        final timestamp = DateTimeFormatter.ISO_INSTANT.format(Instant.now().atOffset(ZoneOffset.UTC))
        sb.append("- Generated: ${timestamp}\n")
        sb.append("- Nextflow version: ${BuildInfo.version}\n")

        // Summary line
        final parts = []
        if( summary.errors > 0 )
            parts.add("${summary.errors} error${summary.errors == 1 ? '' : 's'}")
        if( summary.warnings > 0 )
            parts.add("${summary.warnings} warning${summary.warnings == 1 ? '' : 's'}")
        if( parts.size() > 0 )
            sb.append("- Summary: ${parts.join(', ')}\n")
        else
            sb.append("- Summary: No issues found\n")

        // Sort entries by filename then position
        final sortedErrors = errors.sort { a, b ->
            final cmp = a.filename <=> b.filename
            if( cmp != 0 ) return cmp
            final lineCmp = a.startLine <=> b.startLine
            if( lineCmp != 0 ) return lineCmp
            return a.startColumn <=> b.startColumn
        }

        final sortedWarnings = warnings.sort { a, b ->
            final cmp = a.filename <=> b.filename
            if( cmp != 0 ) return cmp
            final lineCmp = a.startLine <=> b.startLine
            if( lineCmp != 0 ) return lineCmp
            return a.startColumn <=> b.startColumn
        }

        // Errors section
        if( sortedErrors.size() > 0 ) {
            sb.append('\n## :x: Errors\n\n')
            for( final entry : sortedErrors ) {
                sb.append(formatEntry('Error', entry))
            }
        }

        // Warnings section
        if( sortedWarnings.size() > 0 ) {
            sb.append('\n## :warning: Warnings\n\n')
            for( final entry : sortedWarnings ) {
                sb.append(formatEntry('Warning', entry))
            }
        }

        println sb.toString().trim()
    }

    private String formatEntry(String type, LintEntry entry) {
        final sb = new StringBuilder()
        sb.append("- ${type}: `${entry.filename}:${entry.startLine}:${entry.startColumn}`: ${entry.message}\n\n")

        // Add code context
        final lines = entry.source.getSource().getReader().readLines()
        if( entry.startLine > 0 && entry.startLine <= lines.size() ) {
            final line = lines[entry.startLine - 1]
            final startCol = Math.max(0, entry.startColumn - 1)
            final endCol = Math.min(line.length(), Math.max(startCol + 1, entry.endColumn - 1))

            sb.append("    ```nextflow\n")
            sb.append("    ${line}\n")

            // Add caret markers
            final caretCount = Math.max(1, endCol - startCol)
            sb.append("    ${' ' * startCol}${'^' * caretCount}\n")
            sb.append("    ```\n\n")
        }

        return sb.toString()
    }
}
