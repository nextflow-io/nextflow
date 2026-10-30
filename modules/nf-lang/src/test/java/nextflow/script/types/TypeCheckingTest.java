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
import org.codehaus.groovy.ast.ClassHelper;
import org.codehaus.groovy.ast.expr.BinaryExpression;
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
import org.junit.jupiter.params.provider.Arguments;
import org.junit.jupiter.params.provider.MethodSource;
import org.junit.jupiter.params.provider.ValueSource;

import static nextflow.script.types.TypeCheckingUtils.getType;
import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.tuple;
import static org.junit.jupiter.params.provider.Arguments.arguments;

/**
 * Specifies how the type checker checks declarations, assignments
 * and return statements.
 *
 * @see nextflow.script.control.TypeCheckingVisitor
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class TypeCheckingTest {

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
     * Get the expression of the last statement in the entry workflow
     * of a script that type-checks without errors.
     */
    private static Expression parseExpression(String contents) {
        var source = parse(contents);
        assertThat(source.getErrorCollector().getErrors()).isNullOrEmpty();
        var main = ((ScriptNode) source.getAST()).getEntry().main;
        assertThat(main).isInstanceOf(BlockStatement.class);
        var statements = ((BlockStatement) main).getStatements();
        var last = statements.get(statements.size() - 1);
        assertThat(last).isInstanceOf(ExpressionStatement.class);
        return ((ExpressionStatement) last).getExpression();
    }

    private static void assertType(Expression exp, Class<?> type) {
        assertThat(getType(exp))
            .matches(cn -> Types.isEqual(cn, ClassHelper.makeCached(type)), "is " + type.getSimpleName());
    }

    @Nested
    class Declarations {

        @ParameterizedTest
        @ValueSource(strings = {
            "String x = 'a'",
            "def String x = 'a'",
            "def hello(String x) { return x }",
        })
        void reportLegacyTypeAnnotation(String source) {
            assertThat(check(source)).containsExactly("Legacy type annotations are not supported with static typing -- use `name: Type` instead");
        }

        @ParameterizedTest
        @ValueSource(strings = {
            "def x: String = 'a'",
            "def hello(x: String) { return x }",
        })
        void acceptTypeAnnotation(String source) {
            assertThat(getErrors(source)).isEmpty();
        }

        @Test
        void reportParameterDefaultValueWithWrongType() {
            var errors = getErrors("""
                params {
                    input: String = 3.14
                }

                workflow {
                    params
                }
                """);
            assertThat(errors)
                .extracting(SyntaxException::getStartLine, SyntaxException::getStartColumn, SyntaxException::getOriginalMessage)
                .containsExactly(tuple(2, 5, "Parameter `input` with type String cannot be assigned to default value with type Float"));
        }

        @Test
        void resolveParamsAsParamsMapWithDeclaredFields() {
            var exp = parseExpression("""
                params {
                    input: String = 'input'
                }

                workflow {
                    params
                }
                """);
            assertType(exp, ParamsMap.class);
            assertThat(getType(exp).getField("input")).isNotNull();
        }

        static Stream<Arguments> variableAssignments() {
            return Stream.of(
                arguments("def x = 'hello'\n", String.class),
                arguments("def x\nx = 'hello'\n", String.class),
                arguments("x = null\nx = 'hello'\n", String.class),
                arguments("def x = 1\nx += 1\n", Integer.class)
            );
        }

        @ParameterizedTest
        @MethodSource("variableAssignments")
        void inferTypeOfVariableDeclarationOrAssignment(String source, Class<?> type) {
            var exp = (BinaryExpression) parseExpression(source);
            assertType(exp.getLeftExpression(), type);
        }

        @Test
        void reportAssignmentToRecordField() {
            assertThat(check("""
                def person = record(name: 'Alice', age: 42)
                person.age = 43
                """)).containsExactly("Record fields cannot be assigned -- records are immutable");
        }

        @Test
        void reportCompoundAssignmentToRecordField() {
            assertThat(check("""
                record Person {
                    name: String
                    age: Integer
                }

                workflow {
                    def person = record(name: 'Alice', age: 42) as Person
                    person.age += 1
                }
                """)).containsExactly("Record fields cannot be assigned -- records are immutable");
        }

        @Test
        void reportAssigningResultOfVoidFunctionCall() {
            assertThat(check("""
                def foo() {
                    println 'hello'
                }

                workflow {
                    foo()
                    def x = foo()
                }
                """)).containsExactly("Cannot assign the result of a call that does not return a value");
        }

        @Test
        void reportAssigningResultOfProcessWithoutOutputs() {
            assertThat(check("""
                nextflow.enable.types = true

                process foo {
                    input:
                    message: String

                    script:
                    ''
                }

                workflow {
                    def result = foo('hello')
                }
                """)).containsExactly("Cannot assign the result of a call that does not return a value");
        }

        @Test
        void reportAssigningResultOfWorkflowWithoutOutputs() {
            assertThat(check("""
                nextflow.enable.types = true

                workflow foo {
                    take:
                    message: String

                    main:
                    println message
                }

                workflow {
                    def result = foo( 'hello' )
                }
                """)).containsExactly("Cannot assign the result of a call that does not return a value");
        }
    }

    @Nested
    class Assignments {

        static Stream<Arguments> invalidAssignments() {
            return Stream.of(
                arguments("def x: String ; x = 42", "Assignment target with type String cannot be assigned to value with type Integer"),
                arguments("def x: List<String> ; x = [42]", "Assignment target with type List<String> cannot be assigned to value with type List<Integer>")
            );
        }

        @ParameterizedTest
        @MethodSource("invalidAssignments")
        void reportAssignmentWithWrongType(String source, String error) {
            assertThat(check(source)).containsExactly(error);
        }

        @ParameterizedTest
        @ValueSource(strings = {
            "def x: List<String> ; x = []",
            "def ch: Channel<Path> = channel.empty()",
        })
        void acceptEmptyValueAssignedToTypedTarget(String source) {
            assertThat(getErrors(source)).isEmpty();
        }

        @ParameterizedTest
        @ValueSource(strings = {
            "def x = 1 ; x += 1",
            "def x = 1 ; x -= 1",
            "def x = 1 ; x *= 2",
            "def x = 1 ; x <<= 2",
            "def s = 'a' ; s += 'b'",
            "def d = 1.h ; d *= 2.0",
        })
        void acceptCompoundAssignment(String source) {
            assertThat(getErrors(source)).isEmpty();
        }

        static Stream<Arguments> invalidCompoundAssignments() {
            return Stream.of(
                arguments("def s = 'a' ; s *= 'b'", "The `*=` operator is not defined for operands with types String and String"),
                arguments("def x = 1 ; x += 'a'", "The `+=` operator is not defined for operands with types Integer and String")
            );
        }

        @ParameterizedTest
        @MethodSource("invalidCompoundAssignments")
        void reportCompoundAssignmentWithUndefinedOperator(String source, String error) {
            assertThat(check(source)).containsExactly(error);
        }
    }

    @Nested
    class Returns {

        @Test
        void inferReturnTypeOfUntypedFunctionAtCallSite() {
            var source = parse("""
                def make() {
                    return channel.of(1, 2, 3)
                }

                workflow {
                    make()
                }
                """);
            var ast = (ScriptNode) source.getAST();
            var fn = ast.getFunctions().get(0);
            var statements = ((BlockStatement) ast.getEntry().main).getStatements();
            var call = ((ExpressionStatement) statements.get(statements.size() - 1)).getExpression();
            assertThat(Types.getName(getType(call))).isEqualTo("Channel<Integer>");
            // the function itself stays untyped, since code generation reads it
            assertThat(ClassHelper.isDynamicTyped(fn.getReturnType())).isTrue();
        }

        @Test
        void reportReturnValueWithWrongType() {
            var errors = getErrors("""
                def hello() -> String {
                    return 42
                }
                """);
            assertThat(errors)
                .extracting(SyntaxException::getStartLine, SyntaxException::getStartColumn, SyntaxException::getOriginalMessage)
                .containsExactly(tuple(2, 5, "Return value with type Integer does not match the declared return type (String)"));
        }

        @Test
        void acceptReturnValueWithDeclaredType() {
            assertThat(getErrors("""
                def hello() -> String {
                    return 'Hello!'
                }
                """)).isEmpty();
        }

        @Test
        void reportImplicitReturnAgainstStatementItWasInferredFrom() {
            var errors = getErrors("""
                def hello() -> String {
                    println 'Hello!'
                }
                """);
            assertThat(errors)
                .extracting(SyntaxException::getStartLine, SyntaxException::getStartColumn, SyntaxException::getOriginalMessage)
                .containsExactly(tuple(2, 5, "Return value with type void does not match the declared return type (String)"));
        }

        @Test
        void reportMissingReturnAtIfStatement() {
            var errors = getErrors("""
                def hello(x: Integer) -> String {
                    if( x > 1 ) {
                        return 'big'
                    }
                }
                """);
            assertThat(errors)
                .extracting(SyntaxException::getStartLine, SyntaxException::getOriginalMessage)
                .containsExactly(tuple(2, "Missing return statement"));
        }

        @Test
        void reportMissingReturnAtLastStatementOfBranch() {
            var errors = getErrors("""
                def hello(x: Integer) -> String {
                    if( x > 1 ) {
                        def y = 'big'
                    } else {
                        return 'small'
                    }
                }
                """);
            assertThat(errors)
                .extracting(SyntaxException::getStartLine, SyntaxException::getOriginalMessage)
                .containsExactly(tuple(3, "Missing return statement"));
        }

        @Test
        void reportMissingReturnOfEmptyBodyAgainstFunctionOrClosure() {
            var errors = getErrors("""
                def hello(x: Integer) -> String { }

                workflow {
                    channel.of(1).map { x -> }
                }
                """);
            assertThat(errors)
                .extracting(SyntaxException::getStartLine, SyntaxException::getStartColumn, SyntaxException::getOriginalMessage)
                .containsExactly(
                    tuple(1, 1, "Missing return statement"),
                    tuple(4, 23, "Missing return statement"));
        }

        @ParameterizedTest
        @ValueSource(strings = {
            "def f(x: Integer) -> String { if( x > 1 ) { return 'a' } }",
            "def f(x: Integer) -> String { if( x > 1 ) { } else { return 'a' } }",
            "def f(x: Integer) -> String { if( x > 1 ) { def y = 'a' } else { return 'b' } }",
            "def f(x: Integer) -> String { def y = 'a' ; if( x > 1 ) { y = 'b' } else { 'c' } }",
            "def f(x: Integer) -> String { def y = 'a' }",
            "def f(x: Integer) -> String { try { 'a' } catch( e: Exception ) { def y = 'b' } }",
            "channel.of(1).map { x -> if( x > 1 ) { x } }",
            "channel.of(1).map { x -> if( x > 1 ) { x } else { def y = 0 } }",
        })
        void reportMissingReturn(String source) {
            assertThat(check(source)).containsExactly("Missing return statement");
        }

        @ParameterizedTest
        @ValueSource(strings = {
            "def f(x: Integer) -> String { if( x > 1 ) { return 'a' } ; return 'b' }",
            "def f(x: Integer) -> String { if( x > 1 ) { 'a' } else if( x > 0 ) { 'b' } else { 'c' } }",
            "def f(x: Integer) -> String { if( x > 1 ) { return 'a' } else { throw new Exception() } }",
            "def f(x: Integer) -> String { try { 'a' } catch( e: Exception ) { 'b' } }",
            "def f(x: Integer) { if( x > 1 ) { println 'a' } }",
            "def f(x: Integer) { if( x > 1 ) { def y = 'a' } }",
            "channel.of(1).map { x -> if( x > 1 ) { x } else { 0 } }",
        })
        void acceptEveryPathReturnsOrNoReturnTypeDeclared(String source) {
            assertThat(getErrors(source)).isEmpty();
        }

        // a bare or implicit void return has no value to match the declared return type
        @ParameterizedTest
        @ValueSource(strings = { "return", "println('Hello!')" })
        void reportFunctionThatReturnsNoValue(String body) {
            assertThat(check("""
                def hello() -> String {
                    %s
                }
                """.formatted(body))).containsExactly("Return value with type void does not match the declared return type (String)");
        }

        static Stream<Arguments> closuresWithoutValue() {
            return Stream.of(
                arguments("[1, 2].inject(0) { acc, v -> return }", "Return value with type void does not match the declared return type (Integer)"),
                arguments("[1, 2].inject(0) { acc, v -> println(v) }", "Return value with type void does not match the declared return type (Integer)"),
                arguments("['a', 'b'].findAll { s -> return }", "Return value with type void does not match the declared return type (Boolean)"),
                arguments("['a', 'b'].findAll { s -> println(s) }", "Return value with type void does not match the declared return type (Boolean)")
            );
        }

        // a bare or implicit void return has no value to match the declared return type
        @ParameterizedTest
        @MethodSource("closuresWithoutValue")
        void reportClosureThatReturnsNoValue(String source, String error) {
            assertThat(check(source)).containsExactly(error);
        }

        @Test
        void doNotCheckReturnValueOfClosureWithVoidSignature() {
            // each() takes a Consumer, which discards the closure result
            assertThat(getErrors("""
                def total = 0
                [1, 2, 3].each { n -> total = total + n }
                """)).isEmpty();
        }

        @Test
        void doNotCheckReturnValueOfClosureWithBooleanSignature() {
            // findAll() takes a Predicate, which coerces the closure result to a boolean
            assertThat(getErrors("['a', 'b', ''].findAll { s -> s }")).isEmpty();
        }
    }

}
