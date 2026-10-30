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
package nextflow.script.formatter;

import java.util.List;

import nextflow.script.control.ScriptParser;
import nextflow.script.control.ScriptResolveVisitor;
import nextflow.script.control.TypeCheckingVisitor;
import nextflow.script.dsl.Types;
import nextflow.test.TestUtils;

import static org.assertj.core.api.Assertions.assertThat;

/**
 * Helpers to format a script and compare it with the expected output.
 */
final class FormatterTestUtils {

    static final FormattingOptions DEFAULT = new FormattingOptions(4, true);

    private FormatterTestUtils() {}

    /**
     * Format a script, optionally after type checking it. Resolve and type
     * errors are allowed, but syntax errors and missing comments are not.
     */
    static String format(String contents, FormattingOptions options, boolean typeCheck) {
        var parser = new ScriptParser();
        var source = parser.parse("main.nf", contents);
        new ScriptResolveVisitor(source, parser.compiler().compilationUnit(), Types.DEFAULT_SCRIPT_IMPORTS, List.of()).visit();
        if( typeCheck )
            new TypeCheckingVisitor(source).visit();
        assertThat(TestUtils.getSyntaxErrors(source)).isEmpty();
        var formatter = new ScriptFormattingVisitor(source, options);
        formatter.visit();
        assertThat(formatter.getMissingComments()).isEmpty();
        return formatter.toString();
    }

    static String format(String contents, FormattingOptions options) {
        return format(contents, options, false);
    }

    static String format(String contents) {
        return format(contents, DEFAULT, false);
    }

    /**
     * Assert that {@code input} is formatted as {@code expected}, and that
     * {@code expected} is already formatted.
     */
    static void assertFormat(String input, String expected, FormattingOptions options) {
        assertThat(format(input, options)).isEqualTo(expected);
        assertThat(format(expected, options)).as("formatting is idempotent").isEqualTo(expected);
    }

    static void assertFormat(String input, String expected) {
        assertFormat(input, expected, DEFAULT);
    }

    static void assertFormatted(String source) {
        assertThat(format(source)).isEqualTo(source);
    }

}
