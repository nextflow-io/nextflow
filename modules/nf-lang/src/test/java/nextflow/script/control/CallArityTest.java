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
package nextflow.script.control;

import java.util.List;

import nextflow.test.TestUtils;
import org.codehaus.groovy.syntax.SyntaxException;
import org.junit.jupiter.api.Test;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.tuple;

/**
 * Specifies the call arity checks of {@link CallArityVisitor}.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class CallArityTest {

    private static List<SyntaxException> check(String contents) {
        return TestUtils.check(new ScriptParser(), contents);
    }

    @Test
    void reportProcessCallWithWrongNumberOfArguments() {
        var errors = check("""
            process hello {
                script:
                \"""
                echo hello
                \"""
            }

            workflow {
                hello('world')
            }
            """);
        assertThat(errors)
            .extracting(SyntaxException::getStartLine, SyntaxException::getStartColumn, SyntaxException::getOriginalMessage)
            .containsExactly(tuple(9, 5, "Incorrect number of call arguments, expected 0 but received 1"));
    }

    @Test
    void reportWorkflowCallWithWrongNumberOfArguments() {
        var errors = check("""
            workflow hello {
                take:
                x
                y

                main:
                println 'hello'
            }

            workflow {
                hello('world')
            }
            """);
        assertThat(errors)
            .extracting(SyntaxException::getStartLine, SyntaxException::getStartColumn, SyntaxException::getOriginalMessage)
            .containsExactly(tuple(11, 5, "Incorrect number of call arguments, expected 2 but received 1"));
    }

    @Test
    void reportAgentCallWithWrongNumberOfArguments() {
        var errors = check("""
            agent reporter {
                model 'openai/gpt-4o'
                instruction 'i'

                input:
                sample: String

                output:
                stdout()

                prompt:
                \"""
                Report on ${sample}.
                \"""
            }

            workflow {
                reporter('a', 'b')
            }
            """);
        assertThat(errors)
            .extracting(SyntaxException::getStartLine, SyntaxException::getStartColumn, SyntaxException::getOriginalMessage)
            .containsExactly(tuple(18, 5, "Incorrect number of call arguments, expected 1 but received 2"));
    }

    @Test
    void acceptAgentCallMatchingDeclaredInputs() {
        var errors = check("""
            agent reporter {
                model 'openai/gpt-4o'
                instruction 'i'

                input:
                sample: String

                output:
                stdout()

                prompt:
                \"""
                Report on ${sample}.
                \"""
            }

            workflow {
                reporter('a')
            }
            """);
        assertThat(errors).isEmpty();
    }

}
