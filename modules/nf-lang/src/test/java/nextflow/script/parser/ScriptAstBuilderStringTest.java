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
package nextflow.script.parser;

import java.nio.file.Path;
import java.util.List;

import nextflow.script.control.ScriptParser;
import nextflow.test.TestUtils;
import org.codehaus.groovy.syntax.SyntaxException;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import static nextflow.test.TestUtils.assertErrors;
import static nextflow.test.TestUtils.tempFile;
import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.tuple;

/**
 * Specifies the syntax errors reported for string literals: unexpected
 * characters in gstrings and invalid escape sequences.
 *
 * @see ScriptAstBuilder
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class ScriptAstBuilderStringTest {

    private static List<SyntaxException> check(String contents) {
        return TestUtils.check(new ScriptParser(), contents);
    }

    @Nested
    class GStrings {

        @Test
        void bareDollarInTripleQuotedGString() {
            assertErrors(check("""
                \"""
                echo "$/"
                \"""
                """),
                tuple(2, 7, "Unexpected character: '$'"));
        }

        @Test
        void bareDollarInDoubleQuotedGString() {
            assertErrors(check("""
                "abc $/ def"
                """),
                tuple(1, 6, "Unexpected character: '$'"));
        }

        @Test
        void invalidEscapeInGString() {
            assertErrors(check("""
                "abc \\q def"
                """),
                tuple(1, 6, "Invalid escape sequence: '\\q'"));
        }

        @Test
        void escapedDollarIsAllowed() {
            assertThat(check("""
                \"""
                echo "\\$/"
                \"""
                """)).isEmpty();
        }

        @Test
        void interpolationIsAllowed() {
            assertThat(check("""
                \"""
                echo "${'hello'}"
                \"""
                """)).isEmpty();
        }
    }

    @Nested
    class EscapeSequences {

        @Test
        void invalidEscapeReportedAtItsPositionInMultilineString() {
            // not at the start of the string
            assertErrors(check("""
                x = '''
                    gsub(/\\./, "0")
                '''
                """),
                tuple(2, 11, "Invalid escape sequence: '\\.'"));
        }

        @Test
        void invalidEscapeInSingleQuotedString() {
            assertErrors(check("""
                x = 'a\\.b'
                """),
                tuple(1, 7, "Invalid escape sequence: '\\.'"));
        }

        @Test
        void unicodeEscapeRequiresSingleUAndFourHexDigits() {
            assertThat(check("""
                x = '\\uZZZZ'
                y = '\\uu00e9'
                """))
                .extracting(SyntaxException::getOriginalMessage)
                .containsExactly("Invalid escape sequence: '\\u'", "Invalid escape sequence: '\\u'");
        }

        @Test
        void invalidEscapeInGStringFragmentStartingWithSlash() {
            assertErrors(check("""
                x = "/opt\\q ${params.foo}"
                """),
                tuple(1, 10, "Invalid escape sequence: '\\q'"));
        }

        @Test
        void invalidEscapeInGStringFragmentAfterInterpolation() {
            assertErrors(check("""
                x = "out ${params.foo}/res\\q.txt"
                """),
                tuple(1, 27, "Invalid escape sequence: '\\q'"));
        }

        @Test
        void everyInvalidEscapeIsReported() {
            assertThat(check("""
                x = 'a\\.b'
                y = 'c\\,d'
                """))
                .extracting(SyntaxException::getStartLine, SyntaxException::getStartColumn)
                .containsExactly(tuple(1, 7), tuple(2, 7));
        }

        @Test
        void validEscapesAndSlashyStringsAreAllowed() {
            assertThat(check("""
                x = 'a\\n\\t\\\\\\u00e9\\101b'
                y = /a\\.b/
                """)).isEmpty();
        }

        @Test
        void invalidEscapeAfterLineContinuationWithCrlf(@TempDir Path root) {
            // a line continuation followed by CRLF spans only one line
            var main = tempFile(root, "main.nf", "x = '''a\\" + "\r\n" + "b\\.c'''" + "\r\n");

            assertErrors(TestUtils.check(new ScriptParser(), List.of(main)),
                tuple(2, 2, "Invalid escape sequence: '\\.'"));
        }
    }

}
