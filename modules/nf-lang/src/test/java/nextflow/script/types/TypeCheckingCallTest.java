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
import java.util.stream.Stream;

import nextflow.script.ast.FeatureFlagNode;
import nextflow.script.ast.ScriptNode;
import nextflow.script.control.ScriptParser;
import nextflow.script.control.ScriptResolveVisitor;
import nextflow.script.control.TypeCheckingVisitor;
import nextflow.script.dsl.Types;
import nextflow.test.TestUtils;
import org.codehaus.groovy.ast.expr.ConstantExpression;
import org.codehaus.groovy.control.SourceUnit;
import org.codehaus.groovy.control.messages.SyntaxErrorMessage;
import org.codehaus.groovy.syntax.SyntaxException;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.Arguments;
import org.junit.jupiter.params.provider.CsvSource;
import org.junit.jupiter.params.provider.MethodSource;
import org.junit.jupiter.params.provider.ValueSource;

import static org.assertj.core.api.Assertions.assertThat;
import static org.junit.jupiter.params.provider.Arguments.arguments;

/**
 * Specifies how the type checker checks function and method calls.
 *
 * @see nextflow.script.control.TypeCheckingVisitor
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class TypeCheckingCallTest {

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

    @Nested
    class Functions {

        // groovy.json is resolved by name, so it must be on the test runtime classpath
        @ParameterizedTest
        @ValueSource(strings = {
            "String.format('%d', 1)",
            "Integer.parseInt('1')",
            "Math.max(1, 2)",
            "groovy.json.JsonOutput.toJson([1, 2])",
        })
        void resolveStaticMethodCall(String source) {
            assertThat(getErrors(source)).isEmpty();
        }

        @Test
        void reportUnrecognizedStaticMethod() {
            assertThat(check("String.nope()")).containsExactly("Unrecognized method `nope` for type String");
        }

        static Stream<Arguments> invalidFunctionCalls() {
            return Stream.of(
                arguments("env()", "Function `env` expects 1 argument(s) but received 0"),
                arguments("env(42)", "Argument with type Integer is not compatible with parameter of type String"),
                arguments("file()", "Function `file` (with multiple signatures) was called with incorrect number of arguments and/or incorrect argument types"),
                arguments("file(42)", "Function `file` (with multiple signatures) was called with incorrect number of arguments and/or incorrect argument types"),
                arguments("file('input.txt', checkIfExist: true)", "Named param `checkIfExist` is not defined"),
                arguments("file('input.txt', checkIfExists: 'true')", "Named param `checkIfExists` expects a Boolean but received a String"),
                arguments("file('input.csv').splitCsv(header: true, skip: 'x')", "Named param `skip` expects a Integer but received a String")
            );
        }

        @ParameterizedTest
        @MethodSource("invalidFunctionCalls")
        void reportInvalidFunctionCall(String source, String error) {
            assertThat(check(source)).containsExactly(error);
        }

        @ParameterizedTest
        @ValueSource(strings = {
            "env('HELLO')",
            "file('input.txt')",
            "file('input.txt', checkIfExists: true)",
            // splitCsv header accepts a boolean or a list of columns
            "file('input.csv').splitCsv(header: true)",
            "file('input.csv').splitCsv(header: ['id', 'fastq'])",
        })
        void acceptValidFunctionCall(String source) {
            assertThat(getErrors(source)).isEmpty();
        }

        private static final String DEFAULT_PARAMS = """
            def greet(name: String, greeting: String = 'Hello') -> String {
                "${greeting}, ${name}!"
            }

            def pad(prefix: String = '>', s: String) -> String {
                prefix + s
            }

            workflow {
                %s
            }
            """;

        static Stream<Arguments> invalidCallsWithDefaultParams() {
            return Stream.of(
                arguments("greet()", "Function `greet` expects 1 to 2 argument(s) but received 0"),
                arguments("greet('world', 'Hi', 'x')", "Function `greet` expects 1 to 2 argument(s) but received 3"),
                arguments("greet(42)", "Argument with type Integer is not compatible with parameter of type String"),
                arguments("pad(42)", "Argument with type Integer is not compatible with parameter of type String")
            );
        }

        @ParameterizedTest
        @MethodSource("invalidCallsWithDefaultParams")
        void reportInvalidCallWithDefaultParams(String call, String error) {
            assertThat(check(DEFAULT_PARAMS.formatted(call))).containsExactly(error);
        }

        @ParameterizedTest
        @ValueSource(strings = { "greet('world')", "greet('world', 'Hi')", "pad('x')" })
        void acceptCallOmittingDefaultParams(String call) {
            assertThat(getErrors(DEFAULT_PARAMS.formatted(call))).isEmpty();
        }

        @Test
        void reportRecordArgumentMissingField() {
            assertThat(check("""
                def hello(s: Sample) {
                }

                record Sample {
                    id: String
                    fastq: Path
                }

                workflow {
                    hello( record(id: '1') )
                }
                """)).containsExactly("Argument with type Record {\n    id: String\n} is not compatible with parameter of type Sample");
        }

        @Test
        void acceptRecordArgumentWithAllFields() {
            assertThat(getErrors("""
                def hello(s: Sample) {
                }

                record Sample {
                    id: String
                    fastq: Path
                }

                workflow {
                    hello( record(id: '1', fastq: file('1.fastq')) )
                }
                """)).isEmpty();
        }

        private static final String TYPE_ARGUMENTS = """
            nextflow.enable.types = true

            def f(xs: List<String>) {
            }

            def g(n: Names) {
            }

            process P {
                input:
                xs: List<String>

                script:
                "echo ${xs}"
            }

            workflow W {
                take:
                s: Channel<Sample>

                main:
                s.view()
            }

            record Sample {
                id: String
                n: Integer
            }

            record Names {
                xs: List<String>
            }

            workflow {
                %s
            }
            """;

        static Stream<Arguments> callsWithWrongTypeArguments() {
            return Stream.of(
                arguments("f( [1, 2] )", "Argument with type List<Integer> is not compatible with parameter of type List<String>"),
                arguments("W( channel.of(record(id: 'a')) )", "Argument with type Channel<Record {\n    id: String\n}> is not compatible with parameter of type Channel<Sample>"),
                arguments("P( channel.of([1, 2]) )", "Argument with type List<Integer> is not compatible with process input of type List<String>"),
                arguments("g( record(xs: [1, 2]) )", "Argument with type Record {\n    xs: List<Integer>\n} is not compatible with parameter of type Names")
            );
        }

        @ParameterizedTest
        @MethodSource("callsWithWrongTypeArguments")
        void reportCallArgumentWithWrongTypeArguments(String call, String error) {
            assertThat(check(TYPE_ARGUMENTS.formatted(call))).containsExactly(error);
        }

        @ParameterizedTest
        @ValueSource(strings = {
            "f( ['a', 'b'] )",
            "f( [] )",
            "W( channel.of(record(id: 'a', n: 1)) )",
            "P( channel.of(['a']) )",
            "P( channel.empty() )",
            "g( record(xs: ['a']) )",
        })
        void acceptCallArgumentWithMatchingTypeArguments(String call) {
            assertThat(getErrors(TYPE_ARGUMENTS.formatted(call))).isEmpty();
        }

        static Stream<Arguments> conflictingTypeArguments() {
            return Stream.of(
                arguments("channel.of(1, 'a')", "Argument with type String is not compatible with parameter of type Integer"),
                arguments("channel.of('a', 1)", "Argument with type Integer is not compatible with parameter of type String"),
                arguments("channel.of(1, 2, 'a')", "Argument with type String is not compatible with parameter of type Integer"),
                arguments("channel.of([1], ['a'])", "Argument with type List<String> is not compatible with parameter of type List<Integer>"),
                arguments("channel.of(1).reduce(0) { a, b -> 'x' }", "Return value with type String does not match the declared return type (Integer)")
            );
        }

        @ParameterizedTest
        @MethodSource("conflictingTypeArguments")
        void reportCallArgumentsThatInferConflictingTypeArguments(String source, String error) {
            assertThat(check(source)).containsExactly(error);
        }

        @ParameterizedTest
        @ValueSource(strings = { "channel.of(1, 2)", "channel.of([1], [2])" })
        void acceptCallArgumentsThatInferSameTypeArgument(String source) {
            assertThat(getErrors(source)).isEmpty();
        }
    }

    @Nested
    class Methods {

        static Stream<Arguments> invalidMethodCalls() {
            return Stream.of(
                arguments("channel.queue()", "Unrecognized method `queue` for namespace `channel`"),
                arguments("channel.value()", "Function `value` expects 1 argument(s) but received 0"),
                arguments("workflow.outputDir.name()", "Unrecognized method `name` for type Path"),
                arguments("workflow.outputDir.resolve()", "Function `resolve` expects 1 argument(s) but received 0"),
                arguments("record(id: '1', n: 1).subMap('id')", "Argument with type String is not compatible with parameter of type Iterable<String>")
            );
        }

        @ParameterizedTest
        @MethodSource("invalidMethodCalls")
        void reportInvalidMethodCall(String source, String error) {
            assertThat(check(source)).containsExactly(error);
        }

        @ParameterizedTest
        @ValueSource(strings = {
            "channel.of()",
            "channel.of(1, 2, 3)",
            "workflow.outputDir.resolve('hello')",
            "record(id: '1', n: 1).subMap(['id'])",
            "def r: Record = record(id: '1'); r.subMap(['id'])",
            "record Sample { id: String }\ndef f(s: Sample) { s.subMap(['id']) }",
        })
        void acceptValidMethodCall(String source) {
            assertThat(getErrors(source)).isEmpty();
        }

        @ParameterizedTest
        @CsvSource(delimiter = '|', quoteCharacter = '"', textBlock = """
            toSortedList { s -> s.size() } | toSortedList
            toSortedList()                 | toSortedList
            count('a')                     | count
            count()                        | count
            max()                          | max
            min()                          | min
            sum()                          | sum
            groupTuple()                   | groupTuple
            collectFile()                  | collectFile
            """)
        void reportLegacyOperatorAsDiscouragedRatherThanUnrecognized(String call, String operator) {
            assertThat(check("channel.of('a')." + call)).containsExactly("Operator `" + operator + "` is discouraged from use with static typing");
        }
    }

}
