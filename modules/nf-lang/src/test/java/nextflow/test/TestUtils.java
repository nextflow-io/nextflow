/*
 * Copyright 2013-2026, Seqera Labs
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
package nextflow.test;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Comparator;
import java.util.List;

import nextflow.config.control.ConfigParser;
import nextflow.script.control.Compiler;
import nextflow.script.control.PhaseAware;
import nextflow.script.control.Phases;
import nextflow.script.control.ScriptParser;
import nextflow.script.control.SeverityAware;
import org.assertj.core.groups.Tuple;
import org.codehaus.groovy.control.SourceUnit;
import org.codehaus.groovy.control.messages.SyntaxErrorMessage;
import org.codehaus.groovy.syntax.SyntaxException;

import static org.assertj.core.api.Assertions.assertThat;

/**
 * Helpers to parse and analyze sources and collect their diagnostics.
 */
public final class TestUtils {

    private TestUtils() {}

    /**
     * Parse and analyze a script and return its errors.
     */
    public static List<SyntaxException> check(ScriptParser parser, String contents) {
        var source = parser.parse("main.nf", contents);
        parser.analyze();
        return getErrors(source);
    }

    public static List<SyntaxException> check(ScriptParser parser, List<Path> files) {
        for( var path : files )
            parser.parse(path.toFile());
        parser.analyze();
        return getErrors(files, parser.compiler());
    }

    /**
     * Parse and analyze a config file and return its errors.
     */
    public static List<SyntaxException> check(ConfigParser parser, String contents) {
        var source = parser.parse("nextflow.config", contents);
        parser.analyze();
        return getErrors(source);
    }

    public static List<SyntaxException> check(ConfigParser parser, List<Path> files) {
        for( var path : files )
            parser.parse(path.toFile());
        parser.analyze();
        return getErrors(files, parser.compiler());
    }

    public static final Comparator<SyntaxException> ERROR_COMPARATOR = Comparator
        .comparingInt(SyntaxException::getStartLine)
        .thenComparingInt(SyntaxException::getStartColumn);

    /**
     * Get the errors of a source file, sorted by position. Soft errors are
     * excluded, since they are reported as warnings (see CmdLint).
     */
    public static List<SyntaxException> getErrors(SourceUnit source) {
        var errorCollector = source.getErrorCollector();
        if( !errorCollector.hasErrors() )
            return List.of();
        return errorCollector.getErrors().stream()
            .filter(e -> e instanceof SyntaxErrorMessage)
            .map(e -> ((SyntaxErrorMessage) e).getCause())
            .filter(e -> !SeverityAware.isSoftError(e))
            .sorted(ERROR_COMPARATOR)
            .toList();
    }

    public static List<SyntaxException> getErrors(List<Path> files, Compiler compiler) {
        return files.stream()
            .map(file -> compiler.getSources().get(file.toUri()))
            .flatMap(source -> getErrors(source).stream())
            .toList();
    }

    public static List<String> getWarnings(SourceUnit source) {
        var errorCollector = source.getErrorCollector();
        if( !errorCollector.hasWarnings() )
            return List.of();
        return errorCollector.getWarnings().stream()
            .map(w -> w.getMessage())
            .toList();
    }

    /**
     * Get the errors of a source file that were reported during parsing.
     */
    public static List<SyntaxException> getSyntaxErrors(SourceUnit source) {
        return getErrors(source).stream()
            .filter(error -> !(error instanceof PhaseAware pa) || pa.getPhase() == Phases.SYNTAX)
            .toList();
    }

    /**
     * Assert the errors by line, column and message, in order.
     */
    public static void assertErrors(List<SyntaxException> errors, Tuple... expected) {
        assertThat(errors)
            .extracting(SyntaxException::getStartLine, SyntaxException::getStartColumn, SyntaxException::getOriginalMessage)
            .containsExactly(expected);
    }

    /**
     * Create a file (and its parent directories) under {@code parent}.
     */
    public static Path tempFile(Path parent, String name, String contents) {
        try {
            var result = parent.resolve(name);
            Files.createDirectories(result.getParent());
            Files.createFile(result);
            if( contents != null )
                Files.writeString(result, contents);
            return result;
        }
        catch( IOException e ) {
            throw new UncheckedIOException(e);
        }
    }

    public static Path tempFile(Path parent, String name) {
        return tempFile(parent, name, null);
    }

}
