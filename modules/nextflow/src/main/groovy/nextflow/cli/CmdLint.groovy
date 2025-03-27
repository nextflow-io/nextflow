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
import com.beust.jcommander.validators.NoValueValidator
import com.beust.jcommander.IParameterValidator
import com.beust.jcommander.ParameterException
import groovy.io.FileType
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import java.util.regex.Pattern
import java.time.Instant
import java.nio.file.*
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

    @Parameter(
        names = ['-exclude'],
        description = 'File pattern to exclude from linting (can be specified multiple times)'
    )
    List<String> excludePatterns = []

    static class OutputFormatValidator implements IParameterValidator {
        @Override
        void validate(String name, String value) throws ParameterException {
            def allowedValues = ['full', 'extended', 'concise', 'json']
            if (!allowedValues.contains(value)) {
                throw new ParameterException("Output format must be one of $allowedValues (found: $value)")
            }
        }
    }

    @Parameter(
        names = ['-output-format'],
        description = 'Output format for lint results [env: NXF_LINT_FORMAT]. Options: full, extended, concise, json',
        validateWith = OutputFormatValidator
    )
    String outputFormat = System.getenv('NXF_LINT_FORMAT') ?: 'full'

    private ScriptParser scriptParser

    private ConfigParser configParser

    static class LintResults {
        String date
        List<LintError> errors = []
        LintSummary summary = new LintSummary()
    }

    static class LintError {
        String filename
        Integer startLine
        Integer startColumn
        String message
    }

    static class LintSummary {
        Integer numErrors = 0
        Integer numFilesWithErrors = 0
        Integer numFilesWithoutErrors = 0
    }

    private LintResults jsonResults = new LintResults()

    private int numErrors = 0
    private int numFilesWithErrors = 0
    private int numFilesWithoutErrors = 0

    @Override
    String getName() { 'lint' }

    @Override
    void run() {
        if( !args )
            throw new AbortOperationException("Error: No input files specified")

        scriptParser = new ScriptParser()
        configParser = new ConfigParser()

        jsonResults.date = Instant.now().toString()

        final term = ansi()
        if(outputFormat != 'json') {
            term.a("Checking Nextflow code..").newline()
            AnsiConsole.out.print(term)
            AnsiConsole.out.flush()
        }

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

        jsonResults.summary.numErrors = numErrors
        jsonResults.summary.numFilesWithErrors = numFilesWithErrors
        jsonResults.summary.numFilesWithoutErrors = numFilesWithoutErrors

        if(outputFormat == 'json'){
            println groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(jsonResults))
        } else {
            final emojis = [
                "🔍 📋",
                "🕵🏻‍♂️ 🔎",
                "🔬 👨🏻‍💻",
                "📑 ☑️",
                "🧾 ✔️"
            ]
            def rnd = new Random()
            term.cursorUp(1).eraseLine().cursorUp(1).eraseLine()
            if( outputFormat == 'concise' ) // extra newline needed if no code being shown
                term.newline()
            term.bold().a("Nextflow code checks complete! ${emojis[rnd.nextInt(emojis.size())]}").reset().newline()
            if(numFilesWithErrors > 0)
                term.fg(Color.RED).a(" ❌ ${numFilesWithErrors} file${numFilesWithErrors==1 ? '':'s'} had ${numErrors} error${numErrors==1 ? '':'s'}").newline()
            if(numFilesWithoutErrors > 0)
                term.fg(Color.BLUE).a(" ✅ ${numFilesWithoutErrors} file${numFilesWithoutErrors==1 ? '':'s'} had no errors").newline()
            if(numFilesWithErrors == 0 && numFilesWithoutErrors == 0)
                term.a(" No files found to process").newline()
            AnsiConsole.out.print(term)
            AnsiConsole.out.flush()
        }

        // If we found errors, throw an exception so that we get a non-zero exit code
        if (numErrors > 0) {
            throw new AbortOperationException()
        }
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

    void parse(File file) {
        final name = file.getName()
        // Skip excluded files
        if (shouldExcludeFile(file)) {
            log.debug "Skipping excluded file ${file}"
            return
        }

        log.debug "Linting file ${file}"
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
        if(outputFormat == 'json')
            return
        final str = ansi().cursorUp(1).eraseLine().a(Attribute.INTENSITY_FAINT).a("Checking: ${file.getPath().replaceFirst(/^\.\//, '')}").reset().newline().toString()
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

    private Ansi getCodeBlock(SourceUnit source, SyntaxErrorMessage message, Ansi term) {
        def lineStart = message.cause.getStartLine()
        def colStart = message.cause.getStartColumn()
        def lineEnd = message.cause.getEndLine()
        def colEnd = message.cause.getEndColumn()
        def lines
        source.getSource().getReader().withCloseable { reader ->
            lines = reader.readLines()
        }

        int context = outputFormat == 'extended' ? 2 : 0
        int fromLine = Math.max(1, lineStart - context)
        int toLine = Math.min(lines.size(), lineEnd + context)

        // Enforce max 5 lines
        if (toLine - fromLine + 1 > 5) {
            if (lineStart <= 3) {
                toLine = fromLine + 4
            } else if (lineEnd >= lines.size() - 2) {
                fromLine = toLine - 4
            } else {
                fromLine = lineStart - 2
                toLine = lineStart + 2
            }
        }

        for (int i = fromLine; i <= toLine; i++) {
            String fullLine = lines[i - 1]
            int start = (i == lineStart) ? colStart - 1 : 0
            int end = (i == lineEnd) ? colEnd - 1 : fullLine.length()

            // Truncate to max 70 characters if needed
            int maxLen = 70
            int lineLen = fullLine.length()
            int windowStart = 0
            if (lineLen > maxLen) {
                if (start < maxLen - 10) {
                    windowStart = 0
                } else if (end > lineLen - 10) {
                    windowStart = lineLen - maxLen
                } else {
                    windowStart = start - 30
                }
            }

            String line = fullLine.substring(windowStart, Math.min(lineLen, windowStart + maxLen))
            int adjStart = Math.max(0, start - windowStart)
            int adjEnd = Math.max(adjStart + 1, Math.min(end - windowStart, line.length()))

            // Line number
            term.fg(Ansi.Color.BLUE).a(String.format("%3d | ", i)).reset()

            if (i == lineStart) {
                // Print line, with error in red
                term.a(Attribute.INTENSITY_FAINT).a(line.substring(0, adjStart)).reset()
                term.fg(Ansi.Color.RED).a(line.substring(adjStart, adjEnd)).reset()
                term.a(Attribute.INTENSITY_FAINT).a(line.substring(adjEnd)).reset().newline()

                // Print carets underneath pointing to the error
                String marker = ' ' * adjStart
                String carets = '^' * Math.max(1, adjEnd - adjStart)
                term.a("    | ")
                    .fg(Ansi.Color.RED).bold().a(marker + carets).reset().newline()
            } else {
                term.a(Attribute.INTENSITY_FAINT).a(line).reset().newline()
            }
        }

        return term
    }

    private void printErrors(SourceUnit source) {
        final errorMessages = source.getErrorCollector().getErrors()
        def term = ansi().cursorUp(1).eraseLine()
        for( final message : errorMessages ) {
            if( message instanceof SyntaxErrorMessage ) {
                numErrors += 1
                final cause = message.getCause()
                if(outputFormat == 'json') {
                    def error = new LintError(
                        filename: source.getName(),
                        startLine: cause.getStartLine(),
                        startColumn: cause.getStartColumn(),
                        message: cause.getOriginalMessage()
                    )
                    jsonResults.errors << error
                } else {
                    term.bold().a("${source.getName().replaceFirst(/^\.\//, '')}").reset()
                    term.a(":${cause.getStartLine()}:${cause.getStartColumn()}: ")
                    term = highlightString(cause.getOriginalMessage(), term)
                    if( outputFormat != 'concise' ) {
                        term.newline()
                        term = getCodeBlock(source, message, term)
                    }
                    term.newline()
                }
            }
        }
        if(outputFormat != 'json') {
            // Extra newline as next status update will chomp back one
            term.fg(Color.DEFAULT).newline()
            AnsiConsole.out.print(term)
            AnsiConsole.out.flush()
        }
    }
}
