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

import java.util.List;

import nextflow.script.ast.FeatureFlagNode;
import nextflow.script.ast.ScriptNode;
import nextflow.script.control.ScriptParser;
import nextflow.script.control.ScriptResolveVisitor;
import nextflow.script.control.TypeCheckingVisitor;
import nextflow.script.dsl.Types;
import nextflow.test.TestUtils;
import org.codehaus.groovy.ast.expr.ConstantExpression;
import org.codehaus.groovy.ast.expr.Expression;
import org.codehaus.groovy.ast.stmt.BlockStatement;
import org.codehaus.groovy.ast.stmt.ExpressionStatement;
import org.codehaus.groovy.control.SourceUnit;
import org.codehaus.groovy.control.messages.SyntaxErrorMessage;
import org.codehaus.groovy.syntax.SyntaxException;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;

import static nextflow.script.types.TypeCheckingUtils.getType;
import static org.assertj.core.api.Assertions.assertThat;

/**
 * Specifies how the type checker checks workflow, process and agent
 * calls and resolves their outputs.
 *
 * @see nextflow.script.control.TypeCheckingVisitor
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class TypeCheckingDataflowCallTest {

    /**
     * Parse, resolve and type-check a script with static typing enabled.
     */
    private static SourceUnit parse(String contents) {
        var parser = new ScriptParser();
        var source = parser.parse("main.nf", contents);
        var ast = (ScriptNode) source.getAST();
        assertThat(ast).isNotNull();
        ast.addFeatureFlag(new FeatureFlagNode("nextflow.enable.types", new ConstantExpression(true)));
        new ScriptResolveVisitor(source, parser.compiler().compilationUnit(), Types.DEFAULT_SCRIPT_IMPORTS, List.of()).visit();
        new TypeCheckingVisitor(source).visit();
        return source;
    }

    /**
     * Get the errors of a script, sorted by position. Unlike
     * {@link TestUtils#getErrors}, soft errors are included.
     */
    private static List<SyntaxException> getErrors(String contents) {
        var errors = parse(contents).getErrorCollector().getErrors();
        if( errors == null )
            return List.of();
        return errors.stream()
            .filter(e -> e instanceof SyntaxErrorMessage)
            .map(e -> ((SyntaxErrorMessage) e).getCause())
            .sorted(TestUtils.ERROR_COMPARATOR)
            .toList();
    }

    private static List<String> check(String contents) {
        return getErrors(contents).stream()
            .map(SyntaxException::getOriginalMessage)
            .toList();
    }

    /**
     * Get the type name of the last statement in the entry workflow
     * of a script that type-checks without errors.
     */
    private static String typeOf(String contents) {
        var source = parse(contents);
        assertThat(source.getErrorCollector().getErrors()).isNullOrEmpty();
        var main = ((ScriptNode) source.getAST()).getEntry().main;
        assertThat(main).isInstanceOf(BlockStatement.class);
        var statements = ((BlockStatement) main).getStatements();
        var last = statements.get(statements.size() - 1);
        assertThat(last).isInstanceOf(ExpressionStatement.class);
        Expression exp = ((ExpressionStatement) last).getExpression();
        return Types.getName(getType(exp));
    }

    @Nested
    class Workflows {

        @Test
        void reportWorkflowCallWithMissingArgument() {
            assertThat(check("""
                nextflow.enable.types = true

                workflow hello {
                    take:
                    messages: Channel<String>

                    main:
                    messages.view()
                }

                workflow {
                    hello()
                }
                """)).containsExactly("Workflow `hello` expects 1 argument(s) but received 0");
        }

        @Test
        void resolveSingleWorkflowOutputAsValue() {
            assertThat(typeOf("""
                nextflow.enable.types = true

                workflow hello {
                    emit:
                    42
                }

                workflow {
                    hello()
                }
                """)).isEqualTo("Value<Integer>");
        }

        @Test
        void resolveNamedWorkflowOutputAsValue() {
            assertThat(typeOf("""
                nextflow.enable.types = true

                workflow hello {
                    emit:
                    foo = 'hello'
                    bar = 42
                }

                workflow {
                    hello().bar
                }
                """)).isEqualTo("Value<Integer>");
        }
    }

    @Nested
    class Processes {

        @Test
        void acceptProcessCallWithChannelArgument() {
            assertThat(getErrors("""
                nextflow.enable.types = true

                process hello {
                    input:
                    message: String

                    script:
                    ''
                }

                workflow hello_flow {
                    take:
                    messages: Channel<String>

                    main:
                    hello( messages )
                }
                """)).isEmpty();
        }

        @Test
        void reportProcessCallWithMissingArgument() {
            assertThat(check("""
                nextflow.enable.types = true

                process hello {
                    input:
                    message: String

                    script:
                    ''
                }

                workflow {
                    hello()
                }
                """)).containsExactly("Process `hello` expects 1 argument(s) but received 0");
        }

        @Test
        void reportProcessCallWithWrongChannelElementType() {
            assertThat(check("""
                nextflow.enable.types = true

                process hello {
                    input:
                    message: String

                    script:
                    ''
                }

                workflow hello_flow {
                    take:
                    messages: Channel<Integer>

                    main:
                    hello( messages )
                }
                """)).containsExactly("Argument with type Integer is not compatible with process input of type String");
        }

        @Test
        void reportProcessCallWithMultipleChannelArguments() {
            assertThat(check("""
                nextflow.enable.types = true

                process hello {
                    input:
                    message: String
                    target: String

                    script:
                    ''
                }

                workflow hello_flow {
                    take:
                    messages: Channel<String>
                    targets: Channel<String>

                    main:
                    hello( messages, targets )
                }
                """)).containsExactly("Process `hello` was called with multiple channel arguments which can lead to non-deterministic behavior -- make sure that at most one argument is a channel and that all other arguments are dataflow values");
        }

        @Test
        void reportRecordArgumentMissingProcessInputField() {
            assertThat(check("""
                nextflow.enable.types = true

                process hello {
                    input:
                    record(id: String, fastq: Path)

                    exec:
                    println '...'
                }

                workflow {
                    hello( record(id: '1') )
                }
                """)).containsExactly("Argument with type Record {\n    id: String\n} is not compatible with process input of type Record {\n    id: String\n    fastq: Path\n}");
        }

        @Test
        void acceptRecordArgumentWithAllProcessInputFields() {
            assertThat(getErrors("""
                nextflow.enable.types = true

                process hello {
                    input:
                    record(id: String, fastq: Path)

                    exec:
                    println '...'
                }

                workflow {
                    hello( record(id: '1', fastq: file('1.fastq')) )
                }
                """)).isEmpty();
        }

        @Test
        void resolveProcessOutputAsValue() {
            assertThat(typeOf("""
                nextflow.enable.types = true

                process hello {
                    input:
                    target: String

                    output:
                    "Hello, $target!"

                    exec:
                    true
                }

                workflow {
                    hello( 'World' )
                }
                """)).isEqualTo("Value<String>");
        }

        @Test
        void resolveProcessOutputAsValueWhenCalledWithNull() {
            assertThat(typeOf("""
                nextflow.enable.types = true

                process hello {
                    input:
                    target: String?

                    output:
                    "Hello, $target!"

                    exec:
                    true
                }

                workflow {
                    hello( null )
                }
                """)).isEqualTo("Value<String>");
        }

        @Test
        void resolveTupleProcessOutput() {
            assertThat(typeOf("""
                nextflow.enable.types = true

                process hello {
                    input:
                    target: String

                    output:
                    tuple(target, "Hello, $target!")

                    exec:
                    true
                }

                workflow {
                    hello( 'World' )
                }
                """)).isEqualTo("Value<Tuple<String, String>>");
        }

        @Test
        void resolveRecordProcessOutput() {
            assertThat(typeOf("""
                nextflow.enable.types = true

                process hello {
                    input:
                    target: String

                    output:
                    record(target: target, message: "Hello, $target!")

                    exec:
                    true
                }

                workflow {
                    hello( 'World' )
                }
                """)).isEqualTo("Value<Record {\n    target: String\n    message: String\n}>");
        }

        @ParameterizedTest
        @CsvSource(delimiter = '|', quoteCharacter = '"', textBlock = """
            'foo'             | Value<String>
            channel.of('foo') | Channel<String>
            """)
        void resolveProcessOutputWithNullArgument(String arg, String type) {
            assertThat(typeOf("""
                nextflow.enable.types = true

                process hello {
                    input:
                    greeting: String
                    target: String?

                    output:
                    "$greeting, $target!"

                    exec:
                    true
                }

                workflow {
                    hello( %s, null )
                }
                """.formatted(arg))).isEqualTo(type);
        }
    }

    @Nested
    class Agents {

        // an agent call is lifted and unwrapped like a process call
        @Test
        void acceptAgentCallWithMatchingInput() {
            assertThat(getErrors("""
                nextflow.enable.types = true

                agent hello {
                    input:
                    message: String

                    output:
                    stdout()

                    prompt:
                    "say ${message}"
                }

                workflow {
                    hello( channel.of('hi') ).view { v -> v.toUpperCase() }
                }
                """)).isEmpty();
        }

        @Test
        void reportAgentCallWithWrongInputType() {
            assertThat(check("""
                nextflow.enable.types = true

                agent hello {
                    input:
                    message: String

                    output:
                    stdout()

                    prompt:
                    "say ${message}"
                }

                workflow {
                    hello( channel.of(42) )
                }
                """)).containsExactly("Argument with type Integer is not compatible with agent input of type String");
        }
    }

}
