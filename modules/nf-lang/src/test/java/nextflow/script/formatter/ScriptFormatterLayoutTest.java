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

import static nextflow.script.formatter.FormatterTestUtils.DEFAULT;
import static nextflow.script.formatter.FormatterTestUtils.assertFormat;

/**
 * Specifies how {@link ScriptFormattingVisitor} re-indents multi-line strings
 * and aligns trailing comments (Harshil alignment).
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class ScriptFormatterLayoutTest {

    private static final FormattingOptions HARSHIL = new FormattingOptions(4, true, true, false, false);

    @Nested
    class MultiLineStrings {

        @Test
        void shiftStringBodyRightWithStatement() {
            assertFormat(
                """
                process FOO {
                  script:
                  \"""
                  echo hi
                    echo there
                  \"""
                }
                """,
                """
                process FOO {
                    script:
                    \"""
                    echo hi
                      echo there
                    \"""
                }
                """,
                DEFAULT
            );
        }

        @Test
        void shiftStringBodyLeftWithStatement() {
            assertFormat(
                """
                process FOO {
                    script:
                      \"""
                      echo hi
                        echo there
                      \"""
                }
                """,
                """
                process FOO {
                    script:
                    \"""
                    echo hi
                      echo there
                    \"""
                }
                """,
                DEFAULT
            );
        }

        // `beta` has no indentation to give up, so it stays put and the
        // relative indentation of the string body is not preserved
        @Test
        void neverShiftStringLinePastItsContent() {
            assertFormat(
                """
                workflow {
                        x = \"""
                          alpha
                beta
                          \"""
                        println(x)
                }
                """,
                """
                workflow {
                    x = \"""
                      alpha
                beta
                      \"""
                    println(x)
                }
                """,
                DEFAULT
            );
        }

        // the closure body moved because the chain was wrapped, not because
        // its indentation changed, so there is nothing meaningful to shift by
        @Test
        void doNotReindentStringWhoseStatementStartsMidLine() {
            assertFormat(
                """
                workflow {
                  channel.of(1).map { v -> v }.view { \"""
                  value: ${v}
                  \""" }
                }
                """,
                """
                workflow {
                    channel.of(1)
                        .map { v -> v }
                        .view {
                            \"""
                  value: ${v}
                  \"""
                        }
                }
                """,
                DEFAULT
            );
        }
    }

    @Nested
    class HarshilAlignment {

        @Test
        void alignTrailingCommentsInLegacyWorkflow() {
            assertFormat(
                """
                workflow WORKFLOW {
                    take:
                    ch_input // channel: [ val(meta), [ fastq ] ]
                    val_reference //  string: the reference key

                    main:
                    ch_out = ch_input

                    emit:
                    ch_out // channel: out
                    versions = ch_out // channel: versions
                    summary = ch_out.map { v -> v } // channel: summary
                }
                """,
                """
                workflow WORKFLOW {
                    take:
                    ch_input      // channel: [ val(meta), [ fastq ] ]
                    val_reference //  string: the reference key

                    main:
                    ch_out = ch_input

                    emit:
                    ch_out                           // channel: out
                    versions = ch_out                // channel: versions
                    summary  = ch_out.map { v -> v } // channel: summary
                }
                """,
                HARSHIL
            );
        }

        @Test
        void alignTrailingCommentsButNotTypeAnnotationsInTypedWorkflow() {
            assertFormat(
                """
                nextflow.enable.types = true

                workflow WORKFLOW {
                    take:
                    samples: Channel<Tuple<String,Path>> // the samples
                    reference: Path // the reference
                    skip: Boolean

                    main:
                    result = samples

                    emit:
                    result: Channel<Tuple<String,Path>> = result // the result
                    count: Integer = 1 // the count
                }
                """,
                """
                nextflow.enable.types = true

                workflow WORKFLOW {
                    take:
                    samples: Channel<Tuple<String, Path>> // the samples
                    reference: Path                       // the reference
                    skip: Boolean

                    main:
                    result = samples

                    emit:
                    result: Channel<Tuple<String, Path>> = result // the result
                    count: Integer = 1                            // the count
                }
                """,
                HARSHIL
            );
        }

        @Test
        void alignTrailingCommentsInTypedProcess() {
            assertFormat(
                """
                nextflow.enable.types = true

                process hello {
                    input:
                    tuple(id: String, infile: Path) // the sample
                    index: Path // the index
                    opts: Map<String,String>? // extra options

                    output:
                    result: Path = file('output.txt') // the result
                    log = stdout() // the log

                    script:
                    'cat input.txt > output.txt'
                }
                """,
                """
                nextflow.enable.types = true

                process hello {
                    input:
                    tuple(id: String, infile: Path) // the sample
                    index: Path                     // the index
                    opts: Map<String, String>?      // extra options

                    output:
                    result: Path = file('output.txt') // the result
                    log = stdout()                    // the log

                    script:
                    'cat input.txt > output.txt'
                }
                """,
                HARSHIL
            );
        }

        @Test
        void alignTrailingCommentsInParamsBlockAndRecord() {
            assertFormat(
                """
                params {
                    input: String // the input
                    save_intermeds: Boolean = false // save intermediates
                    outdir: String = 'results'
                }

                record FastqPair {
                    id: String // the sample id
                    fastq_1: Path // the first read
                    fastq_2: Path? // the second read
                }

                workflow {
                    println(params.input)
                }
                """,
                """
                params {
                    input: String                   // the input
                    save_intermeds: Boolean = false // save intermediates
                    outdir: String = 'results'
                }

                record FastqPair {
                    id: String     // the sample id
                    fastq_1: Path  // the first read
                    fastq_2: Path? // the second read
                }

                workflow {
                    println(params.input)
                }
                """,
                HARSHIL
            );
        }

        @Test
        void doNotAlignTrailingCommentsOfMultiLineEntries() {
            assertFormat(
                """
                nextflow.enable.types = true

                process hello {
                    input:
                    record(id: String, infile: Path) // the sample
                    index: Path // the index

                    script:
                    'cat input.txt > output.txt'
                }
                """,
                """
                nextflow.enable.types = true

                process hello {
                    input:
                    record(
                        id: String,
                        infile: Path
                    ) // the sample
                    index: Path // the index

                    script:
                    'cat input.txt > output.txt'
                }
                """,
                HARSHIL
            );
        }
    }

}
