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

package nextflow.script.parser.v2

import groovy.transform.CompileStatic
import groovy.transform.Memoized
import org.codehaus.groovy.control.SourceUnit
import org.codehaus.groovy.control.messages.WarningMessage
import org.codehaus.groovy.syntax.SyntaxException
import org.fusesource.jansi.Ansi
import org.fusesource.jansi.AnsiConsole
/**
 * Interfaces and stdout implementation for reporting compilation errors.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
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


class ErrorSummary {
    int errors = 0
    int filesWithErrors = 0
    int filesWithoutErrors = 0
    int filesFormatted = 0
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
        term.reset().fg(Ansi.Color.RED).a("Error").reset()
        term.bold().a(" ${filename}").reset()
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
        term.reset().fg(Ansi.Color.YELLOW).a("Warn").reset()
        term.bold().a("  ${filename}").reset()
        term.a(":${token.getStartLine()}:${token.getStartColumn()}: ")
        term = highlightString(warning.getMessage(), term)
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

            // Left border
            if(i == toLine && i !== startLine) term.fg(color).a("╰").reset().a(" ")
            else term.fg(color).a("│").reset().a(" ")

            // Line number
            term.fg(Ansi.Color.BLUE).a(String.format("%3d | ", i)).reset()

            if( i == startLine ) {
                // Print line with range highlighted
                term.a(Ansi.Attribute.INTENSITY_FAINT).a(line.substring(0, adjStart)).reset()
                term.fg(color).a(line.substring(adjStart, adjEnd)).reset()
                term.a(Ansi.Attribute.INTENSITY_FAINT).a(line.substring(adjEnd)).reset().newline()

                // Left border
                if(i == toLine) term.fg(color).a("╰").reset().a(" ")
                else term.fg(color).a("│").reset().a(" ")

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
        term.cursorUp(1).eraseLine()
        // print extra newline if no code is being shown
        if( mode == 'concise' )
            term.newline()
        term.bold().a("Nextflow linting complete!").reset().newline()
        if( summary.filesWithErrors > 0 ) {
            term.fg(Ansi.Color.RED).a(" ❌ ${summary.filesWithErrors} file${summary.filesWithErrors==1 ? '' : 's'} had ${summary.errors} error${summary.errors==1 ? '' : 's'}").reset().newline()
        }
        if( summary.filesWithoutErrors > 0 ) {
            term.fg(Ansi.Color.GREEN).a(" ✅ ${summary.filesWithoutErrors} file${summary.filesWithoutErrors==1 ? '' : 's'} had no errors")
            if( summary.filesFormatted > 0 )
                term.fg(Ansi.Color.BLUE).a(" (${summary.filesFormatted} formatted)")
            term.reset().newline()
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
