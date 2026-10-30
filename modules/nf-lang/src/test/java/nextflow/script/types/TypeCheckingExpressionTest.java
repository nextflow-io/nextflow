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

import nextflow.script.ast.ASTNodeMarker;
import nextflow.script.control.TypeCheckingVisitor;
import nextflow.script.dsl.Types;
import org.codehaus.groovy.ast.MethodNode;
import org.codehaus.groovy.ast.expr.Expression;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.Arguments;
import org.junit.jupiter.params.provider.MethodSource;
import org.junit.jupiter.params.provider.ValueSource;

import static nextflow.script.types.TypeCheckingTestUtils.assertType;
import static nextflow.script.types.TypeCheckingTestUtils.check;
import static nextflow.script.types.TypeCheckingTestUtils.getErrors;
import static nextflow.script.types.TypeCheckingTestUtils.parseExpression;
import static nextflow.script.types.TypeCheckingTestUtils.typeOf;
import static nextflow.script.types.TypeCheckingUtils.getType;
import static org.assertj.core.api.Assertions.assertThat;
import static org.junit.jupiter.params.provider.Arguments.arguments;

/**
 * Specifies how the type checker checks and resolves operators,
 * literals, closures and property expressions.
 *
 * @see nextflow.script.control.TypeCheckingVisitor
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class TypeCheckingExpressionTest {

    /**
     * Get the type name of the last parameter of the method that
     * a method call resolved to.
     */
    private static String lastParameterType(Expression call) {
        MethodNode method = call.getNodeMetaData(ASTNodeMarker.METHOD_TARGET);
        var params = method.getParameters();
        return Types.getName(params[params.length - 1].getType());
    }

    @Nested
    class Operators {

        static Stream<Arguments> invalidBinaryExpressions() {
            return Stream.of(
                arguments("2 + '2'", "The `+` operator is not defined for operands with types Integer and String"),
                arguments("2 == '2'", "The `==` operator is not defined for operands with types Integer and String"),
                arguments("2 in ['2']", "The `in` operator is not defined for operands with types Integer and List<String>"),
                arguments("[2] == ['2']", "The `==` operator is not defined for operands with types List<Integer> and List<String>"),
                arguments("[[1]] == [['a']]", "The `==` operator is not defined for operands with types List<List<Integer>> and List<List<String>>"),
                arguments("[a: 1] == [a: 'x']", "The `==` operator is not defined for operands with types Map<String, Integer> and Map<String, String>"),
                arguments("[a: 1] == [(1): 1]", "The `==` operator is not defined for operands with types Map<String, Integer> and Map<Integer, Integer>"),
                arguments("tuple(1, 2) == tuple('a', 2)", "The `==` operator is not defined for operands with types Tuple<Integer, Integer> and Tuple<String, Integer>"),
                arguments("tuple(1, 2) == tuple(1, 2, 3)", "The `==` operator is not defined for operands with types Tuple<Integer, Integer> and Tuple<Integer, Integer, Integer>"),
                arguments("record(id: 'a') == [id: 'a']", "The `==` operator is not defined for operands with types Record {\n    id: String\n} and Map<String, String>"),
                arguments("record A { id: String }\ndef f(a: A, m: Map<String,String>) { a == m }", "The `==` operator is not defined for operands with types A and Map<String, String>")
            );
        }

        @ParameterizedTest
        @MethodSource("invalidBinaryExpressions")
        void reportBinaryExpressionWithIncompatibleOperands(String source, String error) {
            assertThat(check(source)).containsExactly(error);
        }

        @ParameterizedTest
        @ValueSource(strings = {
            "2 + 2",
            "2 == 2",
            "2 in [2]",
            "[2] == [2]",
            "[1] == [1.0]",
            "[a: 1] == [a: 1.0]",
            "tuple(1, 2) == tuple(1.0, 2)",
            "record(id: 'a') == record(name: 'b')",
            "record(id: 'a') == record(id: 1)",
            "record Sample { id: String }\ndef f(s: Sample) { s == record(id: 'a') }",
            "record Sample { id: String }\ndef f(s: Sample) { record(id: 'a') != s }",
            "record Sample { id: String }\ndef f(s: Sample) { s == record(name: 'a') }",
            "record Sample { id: String }\ndef f(s: Sample) { s == record(id: 'a', extra: 1) }",
            "record A { id: String }\nrecord B { name: String }\ndef f(a: A, b: B) { a == b }",
            "record A { id: String }\nrecord B { name: String }\ndef f(a: List<A>, b: List<B>) { a == b }",
            "record Sample { id: String }\ndef f(s: Sample) { s == record(id: 1) }",
            "['a', 'b', 'c'][0]",
            "['a', 'b', 'c'][0..-2]",
            "['a', 'b', 'c'][0..1]",
            "['a', 'b', 'c'][1..<3]",
            "file('ab.txt') << 'a' << 'b'",
        })
        void acceptBinaryExpressionWithCompatibleOperands(String source) {
            assertThat(getErrors(source)).isEmpty();
        }

        static Stream<Arguments> binaryExpressionTypes() {
            return Stream.of(
                arguments("1 <=> 2", "Integer"),
                arguments("'a' <=> 'b'", "Integer"),
                arguments("1 < 2", "Boolean"),
                arguments("file('a.txt') << 'a'", "Path"),
                // list index vs slice
                arguments("['a', 'b', 'c'][0]", "String"),
                arguments("['a', 'b', 'c'][0..-2]", "List<String>")
            );
        }

        @ParameterizedTest
        @MethodSource("binaryExpressionTypes")
        void resolveBinaryExpression(String source, String type) {
            assertThat(typeOf(source)).isEqualTo(type);
        }

        @Test
        void reportUnaryExpressionWithUndefinedOperator() {
            assertThat(check("-'42'")).containsExactly("The `-` operator is not defined for an operand with type String");
        }

        @Test
        void acceptUnaryExpression() {
            assertThat(getErrors("-42")).isEmpty();
        }

        static Stream<Arguments> inconsistentTernaries() {
            return Stream.of(
                arguments("true ? 42 : '42'", "Conditional expression has inconsistent types -- true branch has type Integer but false branch has type String"),
                arguments("true ? ['v'] : [:]", "Conditional expression has inconsistent types -- true branch has type List<String> but false branch has type Map<K, V>"),
                arguments("(true ? [1] : 'a').size()", "Conditional expression has inconsistent types -- true branch has type List<Integer> but false branch has type String")
            );
        }

        @ParameterizedTest
        @MethodSource("inconsistentTernaries")
        void reportTernaryWithInconsistentBranchTypes(String source, String error) {
            assertThat(check(source)).containsExactly(error);
        }

        @ParameterizedTest
        @ValueSource(strings = { "true ? 42 : null", "def f(d) { (true ? d : 'a').size() }" })
        void acceptTernaryWithNullOrDynamicBranch(String source) {
            assertThat(getErrors(source)).isEmpty();
        }

        static Stream<Arguments> ternaryTypes() {
            return Stream.of(
                arguments("true ? 42 : null", "Integer?"),
                arguments("true ? [k: 'v'] : [:]", "Map<String, String>"),
                arguments("true ? [:] : [k: 'v']", "Map<String, String>"),
                arguments("true ? ['v'] : []", "List<String>"),
                arguments("true ? [] : ['v']", "List<String>")
            );
        }

        @ParameterizedTest
        @MethodSource("ternaryTypes")
        void resolveTernaryExpression(String source, String type) {
            assertThat(typeOf(source)).isEqualTo(type);
        }

        @Test
        void resolveElvisExpression() {
            var exp = parseExpression("""
                workflow {
                    def meta: Map<String,String> = [:]
                    meta.id ?: 'id'
                }
                """);
            assertType(exp, String.class);
        }

        @Test
        void reportInvalidCast() {
            assertThat(check("42 as List")).containsExactly("Value of type Integer cannot be cast to List");
        }

        @ParameterizedTest
        @ValueSource(strings = { "[] as List<Path>", "'24 h' as Duration", "'8 GB' as MemoryUnit" })
        void acceptValidCast(String source) {
            assertThat(getErrors(source)).isEmpty();
        }

        static Stream<Arguments> invalidSpreadDots() {
            return Stream.of(
                arguments("'*.txt'*.dirName", "Spread-dot is only supported for Iterable types"),
                arguments("files('*.txt')*.dirName", "Unrecognized property `dirName` for element type Path"),
                arguments("'*.txt'*.getDirName()", "Spread-dot is only supported for Iterable types"),
                arguments("files('*.txt')*.getDirName()", "Unrecognized method `getDirName` for element type Path")
            );
        }

        @ParameterizedTest
        @MethodSource("invalidSpreadDots")
        void reportInvalidSpreadDot(String source, String error) {
            assertThat(check(source)).containsExactly(error);
        }

        @ParameterizedTest
        @ValueSource(strings = { "files('*.txt')*.name", "files('*.txt')*.toUriString()" })
        void acceptSpreadDotOverIterable(String source) {
            assertThat(getErrors(source)).isEmpty();
        }
    }

    @Nested
    class Literals {

        @Test
        void recognizeDurationLiterals() {
            assertThat(getErrors("""
                (1.d + 4.h + 15.m + 30.s + 500.ms).toMillis()
                """)).isEmpty();
        }

        @Test
        void recognizeMemoryUnitLiterals() {
            assertThat(getErrors("""
                (1.GB + 500.MB + 4.KB + 64.B).toBytes()
                """)).isEmpty();
        }

        @ParameterizedTest
        @ValueSource(strings = { "[1, 2, '3']", "[1, 'a'].collect { v -> v.size() }" })
        void reportListWithInconsistentElementTypes(String source) {
            assertThat(check(source)).containsExactly("List expression has inconsistent element types -- some elements have type Integer while others have type String");
        }

        @Test
        void acceptListWithConsistentElementTypes() {
            assertThat(getErrors("[1, 2, 3]")).isEmpty();
        }

        static Stream<Arguments> genericTypes() {
            return Stream.of(
                arguments("[1, 2, 3]", "List<Integer>"),
                arguments("[x: 1, y: 2]", "Map<String, Integer>"),
                arguments("[x: '1', y: 2]", "Map<String, ?>"),
                arguments("([x: '1'] + [y: '2']).keySet()", "Set<String>"),
                arguments("([x: '1'] + [y: '2']).values()", "Bag<String>"),
                arguments("[1, 2, 3].first()", "Integer"),
                arguments("channel.value(\"42\")", "Value<String>"),
                arguments("channel.of(1, 2, 3)", "Channel<Integer>"),
                arguments("[1, 2, 3].collect { v -> v * 2 }", "Iterable<Integer>"),
                arguments("channel.of(1).map { v -> v * 2 }", "Channel<Integer>"),
                arguments("channel.fromList(1..10).flatMap { n -> 1..n }", "Channel<Integer>")
            );
        }

        @ParameterizedTest
        @MethodSource("genericTypes")
        void resolveGenericTypes(String source, String type) {
            assertThat(typeOf(source)).isEqualTo(type);
        }
    }

    @Nested
    class Closures {

        @ParameterizedTest
        @ValueSource(strings = {
            "[1, 2, 3].collect { v -> v * 2 }",
            "channel.of(1, 2, 3).map { v -> v * 2 }",
        })
        void resolveGenericFunctionalParameterType(String source) {
            assertThat(lastParameterType(parseExpression(source))).isEqualTo("(Integer) -> Integer");
        }

        @Test
        void resolveToSortedWithTransformClosure() {
            var exp = parseExpression("""
                [1, 2, 3].toSorted { v -> -v }
                """);
            assertThat(lastParameterType(exp)).isEqualTo("(Integer) -> Integer");
            assertThat(Types.getName(getType(exp))).isEqualTo("List<Integer>");
        }

        @Test
        void resolveToSortedWithComparatorClosure() {
            var exp = parseExpression("""
                [1, 2, 3].toSorted { a, b -> a <=> b }
                """);
            assertThat(lastParameterType(exp)).isEqualTo("(Integer, Integer) -> Integer");
            assertThat(Types.getName(getType(exp))).isEqualTo("List<Integer>");
        }

        static Stream<Arguments> comparatorMethods() {
            return Stream.of(
                arguments("toSorted", "List<Tuple<String, Integer>>"),
                arguments("max", "Tuple<String, Integer>"),
                arguments("min", "Tuple<String, Integer>")
            );
        }

        @ParameterizedTest
        @MethodSource("comparatorMethods")
        void preferComparatorVariantOverTupleDestructuring(String method, String type) {
            var source = """
                def pairs = [tuple('b', 2), tuple('a', 3), tuple('c', 1)]
                pairs.%s { a, b -> a[1] <=> b[1] }
                """.formatted(method);
            var exp = parseExpression(source);
            assertThat(lastParameterType(exp)).isEqualTo("(Tuple<String, Integer>, Tuple<String, Integer>) -> Integer");
            assertThat(Types.getName(getType(exp))).isEqualTo(type);
            assertThat(getErrors(source)).isEmpty();
        }

        @Test
        void treatWildcardTypeArgumentAsUnknown() {
            // channel.topic() returns Channel<?>, so the element type is not known
            assertThat(getErrors("channel.topic('versions').view { v -> v.name }")).isEmpty();
        }

        @Test
        void checkClosureParameterOfKnownElementType() {
            assertThat(check("channel.of('a').view { v -> v.name }")).containsExactly("Unrecognized property `name` for type String");
        }
    }

    @Nested
    class Properties {

        static Stream<Arguments> invalidProperties() {
            return Stream.of(
                arguments("workflow.output_dir", "Unrecognized property `output_dir` for namespace `workflow`"),
                arguments("workflow.outputDir.fileName", "Unrecognized property `fileName` for type Path"),
                arguments("[a: 1].entrySet().each { e -> e.key.nope }", "Unrecognized property `nope` for type String"),
                arguments("[a: 1].entrySet().each { e -> e.value.nope }", "Unrecognized property `nope` for type Integer"),
                arguments("[a: 1].entrySet().each { e -> e.nope }", "Unrecognized property `nope` for type Entry<String, Integer>")
            );
        }

        @ParameterizedTest
        @MethodSource("invalidProperties")
        void reportUnrecognizedProperty(String source, String error) {
            assertThat(check(source)).containsExactly(error);
        }

        @ParameterizedTest
        @ValueSource(strings = {
            "workflow.outputDir",
            "workflow.outputDir.name",
            // getter properties with a generic type
            "[a: 1].entrySet().each { e -> def k: String = e.key }",
            "[a: 1].entrySet().each { e -> def v: Integer = e.value }",
        })
        void acceptValidProperty(String source) {
            assertThat(getErrors(source)).isEmpty();
        }

        @Test
        void resolveUntypedMapPropertyAsObject() {
            var exp = parseExpression("""
                workflow {
                    def meta: Map = [:]
                    meta.id
                }
                """);
            assertType(exp, Object.class);
        }

        @Test
        void resolveUntypedMapPropertyAsObjectAfterReassignment() {
            var exp = parseExpression("""
                workflow {
                    def meta: Map = [:]
                    meta = meta + [single_end: true]
                    meta.id
                }
                """);
            assertType(exp, Object.class);
        }
    }

}
