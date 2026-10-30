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
 * Specifies the syntax errors reported when building the AST of a process.
 *
 * @see ScriptAstBuilder
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class ScriptAstBuilderProcessTest {

    private static List<SyntaxException> check(String contents) {
        return TestUtils.check(new ScriptParser(), contents);
    }

    private static void assertErrors(List<SyntaxException> errors, Tuple... expected) {
        assertThat(errors)
            .extracting(SyntaxException::getStartLine, SyntaxException::getStartColumn, SyntaxException::getOriginalMessage)
            .containsExactly(expected);
    }

    @Nested
    class Sections {

        @Test
        void implicitScriptSectionWithOtherSections() {
            assertErrors(check("""
                process hello {
                    input:
                    val greeting

                    \"""
                    echo '${greeting}!'
                    \"""
                }
                """),
                tuple(5, 5, "The `script:` or `exec:` label is required when other sections are present"));
        }

        @Test
        void assignmentAsDirective() {
            assertErrors(check("""
                process hello {
                    container = 'ubuntu'

                    script:
                    ""
                }
                """),
                tuple(2, 5, "Invalid process directive"));
        }

        @Test
        void ifStatementAsDirective() {
            assertErrors(check("""
                process hello {
                    if( true )
                        container = 'ubuntu'
                    else
                        container = 'centos'

                    script:
                    ""
                }
                """),
                tuple(2, 5, "Invalid process directive"));
        }
    }

    @Nested
    class TypedAndLegacyProcesses {

        @Test
        void typedProcessWithoutPreviewFlag() {
            // only the first two errors are asserted
            assertThat(check("""
                process hello {
                    input:
                    message: String

                    output:
                    result: String

                    exec:
                    result = message
                }
                """))
                .extracting(SyntaxException::getStartLine, SyntaxException::getStartColumn, SyntaxException::getOriginalMessage)
                .startsWith(
                    tuple(3, 5, "Typed input declaration is not allowed in legacy process -- set `nextflow.enable.types = true` to use typed processes in this script"),
                    tuple(6, 5, "Typed output declaration is not allowed in legacy process -- set `nextflow.enable.types = true` to use typed processes in this script"));
        }

        @Test
        void typedProcessWithPreviewFlag() {
            assertThat(check("""
                nextflow.enable.types = true

                process hello {
                    input:
                    message: String

                    output:
                    result: String

                    exec:
                    result = message
                }
                """)).isEmpty();
        }

        @Test
        void legacyProcessWithPreviewFlag() {
            // only the first error is asserted
            assertThat(check("""
                nextflow.enable.types = true

                process hello {
                    input:
                    val message

                    output:
                    val result

                    exec:
                    result = message
                }
                """))
                .extracting(SyntaxException::getStartLine, SyntaxException::getStartColumn, SyntaxException::getOriginalMessage)
                .startsWith(tuple(5, 5, "Invalid input declaration in typed process"));
        }

        @Test
        void legacyProcessWithoutPreviewFlag() {
            assertThat(check("""
                process hello {
                    input:
                    val message

                    output:
                    val result

                    exec:
                    result = message
                }
                """)).isEmpty();
        }
    }

    @Nested
    class TypedSections {

        @Test
        void compoundAssignmentAsOutput() {
            assertErrors(check("""
                nextflow.enable.types = true

                process hello {
                    output:
                    target <<= 'world'

                    script:
                    ""
                }
                """),
                tuple(5, 5, "Invalid output declaration in typed process -- must be a name, assignment, or expression"));
        }

        @Test
        void propertyAssignmentAsOutput() {
            assertErrors(check("""
                nextflow.enable.types = true

                process hello {
                    output:
                    target.name = 'world'

                    script:
                    ""
                }
                """),
                tuple(5, 5, "Invalid output declaration in typed process -- must be a name, assignment, or expression"));
        }

        @Test
        void tupleInputWithOneComponent() {
            assertErrors(check("""
                nextflow.enable.types = true

                process hello {
                    input:
                    tuple(id: String)

                    script:
                    ""
                }
                """),
                tuple(5, 5, "Process tuple input must have more than one component"));
        }

        @Test
        void tupleInputWithMultipleComponents() {
            assertThat(check("""
                nextflow.enable.types = true

                process hello {
                    input:
                    tuple(id: String, fastq: Path)

                    script:
                    ""
                }
                """)).isEmpty();
        }

        @Test
        void multilineRecordInputWithTrailingComma() {
            assertThat(check("""
                nextflow.enable.types = true

                process hello {
                    input:
                    record(
                        id: String,
                        fastq: Path,
                    )

                    script:
                    ""
                }
                """)).isEmpty();
        }

        @Test
        void topicAssignment() {
            assertErrors(check("""
                nextflow.enable.types = true

                process hello {
                    topic:
                    versions = stdout()

                    script:
                    ""
                }
                """),
                tuple(5, 5, "Invalid process topic statement"));
        }

        @Test
        void topicStatement() {
            assertThat(check("""
                nextflow.enable.types = true

                process hello {
                    topic:
                    stdout() >> 'versions'

                    script:
                    ""
                }
                """)).isEmpty();
        }
    }

}
