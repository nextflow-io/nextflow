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
package nextflow.config.parser;

import java.util.List;

import nextflow.config.control.ConfigParser;
import nextflow.test.TestUtils;
import org.codehaus.groovy.syntax.SyntaxException;
import org.junit.jupiter.api.Test;

import static nextflow.test.TestUtils.assertErrors;
import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.tuple;

/**
 * Specifies the syntax errors reported by {@link ConfigAstBuilder}.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class ConfigAstBuilderTest {

    private static List<SyntaxException> check(String contents) {
        return TestUtils.check(new ConfigParser(), contents);
    }

    @Test
    void reportInvalidSyntax() {
        var errors = check("""
            process {
                cpus = 4
                memory 8.GB
            }
            """);
        assertErrors(errors, tuple(3, 12, "Unexpected input: '8'"));
    }

    @Test
    void reportIfStatementAtTopLevel() {
        var errors = check("""
            if( true ) {
                process.cpus = 8
            }
            else {
                process.cpus = 4
            }
            """);
        assertErrors(errors, tuple(1, 1, "If statements cannot be mixed with config statements"));
    }

    @Test
    void reportIfStatementInConfigBlock() {
        var errors = check("""
            process {
                withName: 'HELLO' {
                    if( true ) {
                        cpus = 8
                    }
                    else {
                        cpus = 4
                    }
                }
            }
            """);
        assertErrors(errors, tuple(3, 9, "If statements cannot be mixed with config statements"));
    }

    @Test
    void reportTryCatchAroundConfigAssignment() {
        var errors = check("""
            try {
                process.cpus = 4
            }
            catch( Exception e ) {
            }
            """);
        assertErrors(errors, tuple(1, 1, "Try-catch blocks cannot be mixed with config statements"));
    }

    @Test
    void reportTryCatchAroundConfigInclude() {
        var errors = check("""
            try {
                includeConfig 'http://example.com/nextflow.config'
            }
            catch( Exception e ) {
            }
            """);
        assertErrors(errors, tuple(1, 1, "Try-catch blocks cannot be mixed with config statements"));
    }

    @Test
    void reportVariableDeclaration() {
        var errors = check("""
            def trace_timestamp = new Date().format('yyyy-MM-dd_HH-mm-ss')
            """);
        assertErrors(errors, tuple(1, 1, "Variable declarations cannot be mixed with config statements"));
    }

    @Test
    void reportUnexpectedEndOfFile() {
        var errors = check("""
            process {
                withName: 'FOO' {
                    ext.args = '--foo'

                withName: 'BAR' {
                    ext.args = '--bar'
                }
            }
            """);
        assertThat(errors)
            .extracting(SyntaxException::getStartLine, SyntaxException::getOriginalMessage)
            .containsExactly(tuple(9, "Unexpected end of file -- check for an unclosed brace, bracket, parenthesis, or string"));
    }

    @Test
    void reportInvalidEscapeSequenceAtExactPosition() {
        var errors = check("""
            params.x = 'a\\.b'
            """);
        assertErrors(errors, tuple(1, 14, "Invalid escape sequence: '\\.'"));
    }

    @Test
    void reportInvalidEscapeSequenceInGStringFragmentAfterInterpolation() {
        var errors = check("""
            params.x = "${params.y}/res\\q.txt"
            """);
        assertErrors(errors, tuple(1, 28, "Invalid escape sequence: '\\q'"));
    }

    @Test
    void acceptValidEscapesAndSlashyStrings() {
        var errors = check("""
            params.x = 'a\\n\\t\\\\\\u00e9\\101b'
            params.y = /a\\.b/
            """);
        assertThat(errors).isEmpty();
    }

}
