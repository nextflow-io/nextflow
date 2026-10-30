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

import java.util.List;

import nextflow.script.control.ScriptParser;
import nextflow.test.TestUtils;
import org.assertj.core.groups.Tuple;
import org.codehaus.groovy.syntax.SyntaxException;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.tuple;

/**
 * Specifies the syntax errors reported when building the AST of a script:
 * invalid syntax, script declarations, workflows and records.
 *
 * Process sections are specified in {@link ScriptAstBuilderProcessTest},
 * string literals in {@link ScriptAstBuilderStringTest}.
 *
 * @see ScriptAstBuilder
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class ScriptAstBuilderTest {

    private static final String UNEXPECTED_EOF = "Unexpected end of file -- check for an unclosed brace, bracket, parenthesis, or string";

    private static List<SyntaxException> check(String contents) {
        return TestUtils.check(new ScriptParser(), contents);
    }

    private static void assertErrors(List<SyntaxException> errors, Tuple... expected) {
        assertThat(errors)
            .extracting(SyntaxException::getStartLine, SyntaxException::getStartColumn, SyntaxException::getOriginalMessage)
            .containsExactly(expected);
    }

    private static void assertErrorRanges(List<SyntaxException> errors, Tuple... expected) {
        assertThat(errors)
            .extracting(
                SyntaxException::getStartLine, SyntaxException::getStartColumn,
                SyntaxException::getEndLine, SyntaxException::getEndColumn,
                SyntaxException::getOriginalMessage)
            .containsExactly(expected);
    }

    @Nested
    class InvalidSyntax {

        @Test
        void missingCommaInArguments() {
            assertErrorRanges(check("""
                workflow {
                    foo(1, 2 3, 4)
                }
                """),
                tuple(2, 14, 2, 15, "Unexpected input: '3'"));
        }

        @Test
        void colonAfterDirective() {
            assertErrorRanges(check("""
                process hello {
                    debug true:

                    script:
                    \"""
                    echo hello
                    \"""
                }
                """),
                tuple(2, 16, 2, 17, "Unexpected input: '\\n'"));
        }

        @Test
        void missingCommaInTupleOutput() {
            assertErrorRanges(check("""
                process hello {
                    output:
                    tuple val("hello") val("goodbye")

                    script:
                    \"""
                    echo hello
                    \"""
                }
                """),
                tuple(3, 24, 3, 27, "Unexpected input: 'val'"));
        }

        @Test
        void missingCommaBeforeEmitOption() {
            assertErrorRanges(check("""
                process hello {
                    output:
                    tuple val('hello'), val('goodbye') emit: message

                    script:
                    \"""
                    echo hello
                    \"""
                }
                """),
                tuple(3, 40, 3, 44, "Unexpected input: 'emit'"));
        }

        @Test
        void trailingCommaInWorkflowTake() {
            assertErrorRanges(check("""
                workflow foo {
                    take:
                    x
                    y,
                    z

                    main:
                    println 'hello'
                }
                """),
                tuple(4, 6, 4, 7, "Unexpected input: ','"));
        }

        @Test
        void unclosedBraceReportedAtEndOfFile() {
            assertThat(check("""
                workflow {
                    if( true ) {
                        println 'hello'
                }
                """))
                .extracting(SyntaxException::getStartLine, SyntaxException::getOriginalMessage)
                .containsExactly(tuple(5, UNEXPECTED_EOF));
        }

        @Test
        void unclosedBracketAtEndOfFile() {
            assertThat(check("""
                workflow {
                    println 'hello'
                }
                x = [1, 2
                """))
                .extracting(SyntaxException::getOriginalMessage)
                .containsExactly(UNEXPECTED_EOF);
        }

        @Test
        void unclosedStringAtEndOfFile() {
            assertThat(check("""
                workflow {
                    println \"""hello
                }
                """))
                .extracting(SyntaxException::getOriginalMessage)
                .containsExactly(UNEXPECTED_EOF);
        }

        @Test
        void quotedProcessName() {
            assertErrors(check("""
                process 'MY-PROC' {
                    script:
                    "echo hello"
                }
                """),
                tuple(1, 9, "Unexpected input: ''MY-PROC''"));
        }

        @Test
        void typedVariableDeclarationInLegacyWorkflow() {
            assertErrors(check("""
                workflow {
                    x = 1
                    y: Integer = 2
                }
                """),
                tuple(3, 6, "Unexpected input: ':'"));
        }

        @Test
        void processInputAfterOutput() {
            assertThat(check("""
                process hello {
                    output:
                    stdout

                    input:
                    val x

                    script:
                    "echo $x"
                }
                """))
                .extracting(SyntaxException::getStartLine)
                .containsExactly(5);
        }
    }

    @Nested
    class ScriptDeclarations {

        @Test
        void statementsMixedWithDeclarations() {
            assertThat(check("""
                println 'Hello world!'

                workflow {
                }
                """))
                .singleElement()
                .satisfies(error -> {
                    assertThat(error.getStartLine()).isEqualTo(1);
                    assertThat(error.getStartColumn()).isEqualTo(1);
                    assertThat(error.getOriginalMessage()).contains("Statements cannot be mixed with script declarations");
                });
        }

        @Test
        void paramsBlockWithoutEntryWorkflow() {
            assertErrors(check("""
                params {
                    greeting: String
                }
                """),
                tuple(1, 1, "Params block cannot be defined without an entry workflow"));
        }

        @Test
        void paramsBlockWithEntryWorkflow() {
            assertThat(check("""
                params {
                    greeting: String
                }

                workflow {
                }
                """)).isEmpty();
        }

        @Test
        void outputBlockWithoutEntryWorkflow() {
            assertErrors(check("""
                output {
                }
                """),
                tuple(1, 1, "Output block cannot be defined without an entry workflow"));
        }

        @Test
        void outputBlockWithEntryWorkflow() {
            assertThat(check("""
                workflow {
                }

                output {
                }
                """)).isEmpty();
        }
    }

    @Nested
    class Workflows {

        @Test
        void unnamedEmitAmongMultipleEmits() {
            assertErrors(check("""
                workflow hello {
                    emit:
                    'hello'
                    target = 'world'
                }
                """),
                tuple(3, 5, "Every emit must be assigned to a name when there are multiple emits"));
        }

        @Test
        void ifStatementAsEmit() {
            assertErrors(check("""
                workflow hello {
                    emit:
                    if( true )
                        'hello'
                    else
                        'world'
                }
                """),
                tuple(3, 5, "Invalid workflow emit -- must be a name, assignment, or expression"));
        }

        @Test
        void compoundAssignmentAsEmit() {
            assertErrors(check("""
                workflow hello {
                    emit:
                    target <<= 'world'
                }
                """),
                tuple(3, 5, "Invalid workflow emit -- must be a name, assignment, or expression"));
        }

        @Test
        void publishStatementWithoutAssignment() {
            assertErrors(check("""
                workflow {
                    main:
                    a = channel.of(1)

                    publish:
                    a
                }

                output {
                    a {
                        path 'a'
                    }
                }
                """),
                tuple(6, 5, "Invalid workflow publish statement -- must be an assignment"));
        }
    }

    @Nested
    class Records {

        @Test
        void missingRecordBody() {
            assertErrors(check("""
                record Sample {}
                """),
                tuple(1, 1, "Missing record body"));
        }

        @Test
        void missingFieldType() {
            assertErrors(check("""
                record Sample {
                    id
                }
                """),
                tuple(2, 5, "Missing field type"));
        }
    }

}
