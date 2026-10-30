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
package nextflow.script.types;

import nextflow.script.control.TypeCheckingVisitor;
import org.codehaus.groovy.syntax.SyntaxException;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;

import static nextflow.script.types.TypeCheckingTestUtils.check;
import static nextflow.script.types.TypeCheckingTestUtils.getErrors;
import static nextflow.test.TestUtils.assertErrors;
import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.tuple;

/**
 * Specifies how the type checker checks workflow emits, process outputs
 * and the output block.
 *
 * @see nextflow.script.control.TypeCheckingVisitor
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class TypeCheckingOutputTest {

    @Nested
    class WorkflowEmits {

        @Test
        void allowSingleNamedEmit() {
            assertThat(getErrors("""
                nextflow.enable.types = true

                workflow hello {
                    emit:
                    result = 42
                }
                """)).isEmpty();
        }

        @Test
        void reportWorkflowEmitWithWrongType() {
            var errors = getErrors("""
                nextflow.enable.types = true

                workflow hello {
                    emit:
                    a: Value<String> = channel.value(42)
                    b: Value<Integer> = channel.value(1)
                }
                """);
            assertErrors(errors, tuple(5, 5, "Assignment target with type Value<String> cannot be assigned to value with type Value<Integer>"));
        }

        @Test
        void acceptWorkflowEmitsWithMatchingTypes() {
            assertThat(getErrors("""
                nextflow.enable.types = true

                workflow hello {
                    emit:
                    a: Value<Integer> = channel.value(42)
                    b: Channel<Integer> = channel.of(1)
                    c: Value<Integer>? = null
                    d = 1
                }
                """)).isEmpty();
        }

        @Test
        void reportWorkflowEmitThatIsNotDataflowType() {
            var errors = getErrors("""
                nextflow.enable.types = true

                workflow hello {
                    emit:
                    a: Integer = 42
                    b: Value<Integer> = channel.value(1)
                }
                """);
            assertErrors(errors, tuple(5, 5, "Workflow emit `a` must be declared as a Channel or Value, not Integer"));
        }

        @Test
        void reportTypedEmitThatDoesNotMatchMainVariable() {
            assertThat(check("""
                nextflow.enable.types = true

                workflow hello {
                    main:
                    x = channel.of('a')

                    emit:
                    x: Channel<Integer>
                }
                """)).containsExactly("Workflow emit `x` with type Channel<Integer> cannot be assigned to value with type Channel<String>");
        }

        @Test
        void reportTypedEmitThatDoesNotMatchTakeVariable() {
            assertThat(check("""
                nextflow.enable.types = true

                workflow hello {
                    take:
                    x: Channel<String>

                    emit:
                    x: Channel<Integer>
                }
                """)).containsExactly("Workflow emit `x` with type Channel<Integer> cannot be assigned to value with type Channel<String>");
        }

        @Test
        void acceptTypedEmitThatMatchesMainVariable() {
            assertThat(getErrors("""
                nextflow.enable.types = true

                workflow hello {
                    main:
                    x = channel.of(1)

                    emit:
                    x: Channel<Integer>
                }
                """)).isEmpty();
        }
    }

    @Nested
    class ProcessOutputs {

        @Test
        void allowSingleNamedOutput() {
            assertThat(getErrors("""
                nextflow.enable.types = true

                process bye {
                    output:
                    sample: String = 'hi'

                    script:
                    ''
                }
                """)).isEmpty();
        }

        @Test
        void reportTypedProcessOutputThatDoesNotMatchBodyVariable() {
            assertThat(check("""
                nextflow.enable.types = true

                process hello {
                    output:
                    x: Integer

                    exec:
                    x = 'a'
                }
                """)).containsExactly("Process output `x` with type Integer cannot be assigned to value with type String");
        }
    }

    @Nested
    class OutputBlock {

        @Test
        void acceptParameterReferenceInOutputBlock() {
            assertThat(getErrors("""
                params {
                    outdir: String = 'results'
                }

                workflow {
                    main:
                    ch = channel.of('a')

                    publish:
                    samples = ch
                }

                output {
                    samples: Channel<String> {
                        path params.outdir
                    }
                }
                """)).isEmpty();
        }

        @Test
        void reportInvalidPublishStatementsInOutputPathDirective() {
            var errors = getErrors("""
                nextflow.enable.types = true

                workflow {
                    main:
                    ch_samples = channel.empty()

                    publish:
                    samples = ch_samples
                }

                output {
                    samples: Channel<Sample> {
                        path { s ->
                            s.id >> 42
                            s.fastq >> params.save_fastq ? 'fastq' : null
                        }
                    }
                }

                record Sample {
                    id: String
                    fastq: Path
                }
                """);
            assertThat(errors).satisfiesExactly(
                e -> assertThat(e)
                    .extracting(SyntaxException::getStartLine, SyntaxException::getStartColumn, SyntaxException::getOriginalMessage)
                    .containsExactly(14, 13, "Publish source should be a Path or Iterable<Path> but was specified as a String"),
                e -> assertThat(e)
                    .extracting(SyntaxException::getStartLine, SyntaxException::getStartColumn, SyntaxException::getOriginalMessage)
                    .containsExactly(14, 21, "Publish target should be a String but was specified as a Integer"),
                e -> {
                    assertThat(e).extracting(SyntaxException::getStartLine, SyntaxException::getStartColumn).containsExactly(15, 13);
                    assertThat(e.getOriginalMessage()).contains("Statement is not a valid publish statement");
                });
        }

        @Test
        void acceptParenthesizedTernaryAsPublishTarget() {
            assertThat(getErrors("""
                nextflow.enable.types = true

                workflow {
                    main:
                    ch_samples = channel.empty()

                    publish:
                    samples = ch_samples
                }

                output {
                    samples: Channel<Sample> {
                        path { s ->
                            s.fastq >> (params.save_fastq ? 'fastq' : null)
                        }
                    }
                }

                record Sample {
                    id: String
                    fastq: Path
                }
                """)).isEmpty();
        }
    }

}
