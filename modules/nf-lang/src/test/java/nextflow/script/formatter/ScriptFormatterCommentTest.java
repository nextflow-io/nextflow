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

import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;

import static nextflow.script.formatter.FormatterTestUtils.assertFormat;
import static nextflow.script.formatter.FormatterTestUtils.assertFormatted;
import static nextflow.script.formatter.FormatterTestUtils.format;
import static org.assertj.core.api.Assertions.assertThat;

/**
 * Specifies how {@link ScriptFormattingVisitor} preserves comments.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class ScriptFormatterCommentTest {

    @Nested
    class TrailingComments {

        @Test
        void collapseSpaceBeforeTrailingComment() {
            assertFormat(
                """
                params.index = null     // path to the index
                """,
                """
                params.index = null // path to the index
                """
            );
        }

        @Test
        void preserveTrailingCommentsInBlock() {
            assertFormatted("""
                workflow {
                    x = 1 // about x
                    y = 2 // about y
                }
                """);
        }

        @Test
        void preserveTrailingCommentOnLastStatement() {
            assertFormatted("""
                workflow {
                    x = 1 // the last statement
                }
                """);
        }

        // the branch is expanded onto its own line, so the comment trails
        // the statement as a whole
        @Test
        void preserveTrailingCommentOnSingleLineIf() {
            assertFormat(
                """
                workflow {
                    if (!group) group = [1] // if new, create a new entry
                }
                """,
                """
                workflow {
                    if (!group) {
                        group = [1]
                    } // if new, create a new entry
                }
                """
            );
        }

        @Test
        void moveCommentAfterClosureParametersIntoBody() {
            assertFormat(
                """
                workflow {
                    x = [1, 2, 2, 3].reduce([:]) { acc, v -> // 'acc' collects all values
                        acc[v] = 1
                        return acc // give it back
                    }
                }
                """,
                """
                workflow {
                    x = [1, 2, 2, 3].reduce([:]) { acc, v ->
                        // 'acc' collects all values
                        acc[v] = 1
                        return acc // give it back
                    }
                }
                """
            );
        }

        @Test
        void preserveCommentsInCrlfSource() {
            var source = "// lead\r\nworkflow {\r\n    x = 1   // trail\r\n}\r\n";
            assertThat(format(source)).isEqualTo("// lead\nworkflow {\n    x = 1 // trail\n}\n");
        }

        @Test
        void hoistCommentThatHasNoPlaceToGo() {
            assertFormat(
                """
                workflow {
                    x = 1 +
                        // about the second operand
                        2
                }
                """,
                """
                workflow {
                    // about the second operand
                    x = 1 + 2
                }
                """
            );
        }

        @Test
        void doNotTreatCommentLikeStringTextAsComment() {
            assertFormatted("""
                workflow {
                    ch = channel.fromFilePairs("${params.dir}/*.vcf.gz")
                    url = "${params.host}//path"
                }
                """);
        }
    }

    @Nested
    class DanglingComments {

        @Test
        void preserveDanglingCommentAtEndOfBlock() {
            assertFormatted("""
                workflow {
                    x = channel.of(1)

                    // PROCESS_1(x)
                    // PROCESS_2(x)
                }
                """);
        }

        @Test
        void preserveDanglingCommentInEmptyBlock() {
            assertFormatted("""
                workflow {
                    // nothing here yet
                }
                """);
        }

        @Test
        void preserveDanglingCommentAtEndOfFile() {
            assertFormatted("""
                workflow {
                    println('hi')
                }

                // process something {
                // }
                """);
        }
    }

    @Nested
    class Declarations {

        @Test
        void preserveBlockComment() {
            assertFormatted("""
                /*
                 * A block comment
                 */
                workflow {
                    println('hi')
                }
                """);
        }

        @Test
        void preserveFileHeader() {
            assertFormatted("""
                #!/usr/bin/env nextflow
                /*
                 * Copyright 2013-2024, Seqera Labs
                 */

                // enable the output DSL
                nextflow.preview.output = true

                /*
                 * Command line input parameter
                 */
                params.query = 'data/sample.fa'
                """);
        }

        @Test
        void preserveCommentsInParamsBlock() {
            assertFormatted("""
                params {
                    // Comma-separated list of IDs.
                    input: String

                    // Whether to save intermediate outputs.
                    save_intermeds: Boolean = false
                }

                workflow {
                    println(params.input)
                }
                """);
        }

        @Test
        void preserveCommentsOnFunction() {
            assertFormatted("""
                /*
                 * does a thing
                 */
                def foo(x) {
                    return x // give it back
                }
                """);
        }

        @Test
        void preserveCommentsInProcess() {
            assertFormatted("""
                process FOO {
                    // a directive
                    cpus 2

                    input:
                    // about x
                    val x

                    output:
                    path 'out.txt' // the output

                    script:
                    \"""
                    echo hi
                    \"""
                }
                """);
        }

        @Test
        void preserveCommentsInAgent() {
            assertFormatted("""
                nextflow.enable.types = true

                agent reporter {
                    // a directive
                    model 'openai/gpt-4o'

                    input:
                    // about sample
                    sample: String

                    output:
                    report: String // the report

                    prompt:
                    \"""
                    Write a report for ${sample}.
                    \"""
                }
                """);
        }
    }

    @Nested
    class Expressions {

        @Test
        void preserveCommentsInListAndMapLiterals() {
            assertFormatted("""
                workflow {
                    x = [
                        // first
                        1,
                        2, // second
                    ]
                    y = [
                        a: 1, // one
                        // two
                        b: 2,
                    ]
                }
                """);
        }

        @Test
        void preserveCommentsInCallArguments() {
            assertFormatted("""
                workflow {
                    FOO(
                        // first arg
                        1,
                        2, // second arg
                    )
                }
                """);
        }

        @Test
        void preserveCommentsInIfElse() {
            assertFormatted("""
                workflow {
                    if (true) {
                        // then branch
                        x = 1
                    }
                    else {
                        // else branch
                        x = 2
                    }
                }
                """);
        }

        @Test
        void preserveCommentBeforeElseBranch() {
            assertFormatted("""
                workflow {
                    if (x) {
                        a()
                    }
                    // otherwise
                    else if (y) {
                        b()
                    }
                    // last resort
                    else {
                        c()
                    }
                }
                """);
        }
    }

    @Nested
    class MethodChains {

        @Test
        void endOfLineCommentTrailsPrecedingLink() {
            assertFormatted("""
                workflow {
                    result = channel.of(1, 2)
                        .map { it * 2 } // double
                        .filter { it > 2 }
                        .view()
                }
                """);
        }

        @Test
        void ownLineCommentLeadsFollowingLink() {
            assertFormatted("""
                workflow {
                    result = channel.of(1, 2)
                        .map { it * 2 }
                        // about the filter
                        .filter { it > 2 }
                        .view()
                }
                """);
        }

        @Test
        void commentForcesChainToWrapAtEachLink() {
            assertFormat(
                """
                workflow {
                    result = channel.of(1, 2).map { it * 2 } // double
                        .view()
                }
                """,
                """
                workflow {
                    result = channel.of(1, 2)
                        .map { it * 2 } // double
                        .view()
                }
                """
            );
        }

        @Test
        void breakNamespaceReceiverOntoOwnLineForComment() {
            assertFormatted("""
                workflow {
                    channel // about of
                        .of(1)
                        .view()
                }
                """);
        }
    }

}
