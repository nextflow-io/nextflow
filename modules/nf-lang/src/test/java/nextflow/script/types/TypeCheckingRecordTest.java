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
import nextflow.script.control.SeverityAware;
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
import org.junit.jupiter.params.provider.Arguments;
import org.junit.jupiter.params.provider.MethodSource;
import org.junit.jupiter.params.provider.ValueSource;

import static nextflow.script.types.TypeCheckingUtils.getType;
import static org.assertj.core.api.Assertions.assertThat;
import static org.junit.jupiter.params.provider.Arguments.arguments;

/**
 * Specifies how the type checker checks and resolves records, tuples,
 * empty collections and channel operators.
 *
 * @see nextflow.script.control.TypeCheckingVisitor
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class TypeCheckingRecordTest {

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

    private static String typeOf(String contents) {
        return Types.getName(getType(parseExpression(contents)));
    }

    @Nested
    class Records {

        @Test
        void resolveRecordType() {
            var type = getType(parseExpression("""
                record(id: '1', count: 42)
                """));
            assertThat(Types.getName(type)).isEqualTo("Record {\n    id: String\n    count: Integer\n}");
            assertThat(type.getField("id")).isNotNull();
            assertThat(type.getField("count")).isNotNull();
        }

        @Test
        void resolveRecordField() {
            assertThat(typeOf("""
                record(id: '1', count: 42).id
                """)).isEqualTo("String");
        }

        @Test
        void resolveRecordSum() {
            var type = getType(parseExpression("""
                record(id: '1') + record(count: 42)
                """));
            assertThat(Types.getName(type)).isEqualTo("Record {\n    id: String\n    count: Integer\n}");
            assertThat(type.getField("id")).isNotNull();
            assertThat(type.getField("count")).isNotNull();
        }

        @Test
        void reportRecordTypeUsedAsConstructor() {
            assertThat(check("""
                record Sample {
                    id: String
                    reads: List<Path>
                }

                workflow {
                    new Sample()
                }
                """)).containsExactly("Record type Sample cannot be used as a constructor -- use `record()` instead");
        }

        @Test
        void acceptRecordCastWithAllFields() {
            assertThat(getErrors("""
                workflow {
                    record(id: '1', fastq_1: file('1_1.fastq'), fastq_2: file('1_2.fastq')) as FastqPair
                }

                record FastqPair {
                    id: String
                    fastq_1: Path
                    fastq_2: Path
                }
                """)).isEmpty();
        }

        @Test
        void reportRecordCastWithMissingField() {
            assertThat(check("""
                workflow {
                    record(id: '1', fastq_1: file('1_1.fastq')) as FastqPair
                }

                record FastqPair {
                    id: String
                    fastq_1: Path
                    fastq_2: Path
                }
                """)).containsExactly("Record mismatch -- source record is missing field `fastq_2` required by FastqPair");
        }

        @Test
        void reportRecordCastWithWrongFieldType() {
            assertThat(check("""
                workflow {
                    record(id: '1', fastq_1: file('1_1.fastq'), fastq_2: '1_2.fastq') as FastqPair
                }

                record FastqPair {
                    id: String
                    fastq_1: Path
                    fastq_2: Path
                }
                """)).containsExactly("Record mismatch -- field `fastq_2` of FastqPair expects a Path but received a String");
        }

        @Test
        void acceptRecordChannelAssignedToRecordTypeChannel() {
            assertThat(getErrors("""
                workflow {
                    def ch: Channel<Sample> = channel.of( record(id: '1', fastq: '1.fastq') )
                }

                record Sample {
                    id: String
                    fastq: String
                }
                """)).isEmpty();
        }

        @Test
        void reportRecordChannelMissingFieldOfRecordTypeChannel() {
            assertThat(check("""
                workflow {
                    def ch: Channel<Sample> = channel.of( record(id: '1') )
                }

                record Sample {
                    id: String
                    fastq: String
                }
                """)).containsExactly("Assignment target with type Channel<Sample> cannot be assigned to value with type Channel<Record {\n    id: String\n}>");
        }
    }

    @Nested
    class Tuples {

        @Test
        void resolveTupleType() {
            assertThat(typeOf("""
                tuple('hello', 42, true)
                """)).isEqualTo("Tuple<String, Integer, Boolean>");
        }

        @Test
        void resolveTupleElement() {
            assertThat(typeOf("""
                tuple('hello', 42, true)[1]
                """)).isEqualTo("Integer");
        }

        @Test
        void resolveTupleAssignment() {
            assertThat(typeOf("""
                (x, y, z) = tuple('hello', 42, true)
                tuple(x, y, z)
                """)).isEqualTo("Tuple<String, Integer, Boolean>");
        }

        @Test
        void allowTupleDestructuringInClosure() {
            assertThat(getErrors("""
                channel.value( tuple(1, 2, 3) ).map { x, y, z -> y }
                """)).isEmpty();
        }
    }

    @Nested
    class EmptyCollections {

        static Stream<Arguments> emptyChannels() {
            return Stream.of(
                arguments("def ch = flag ? channel.of(record(id: '1')) : channel.empty()", "missing", "Unrecognized property `missing` for type Record {\n    id: String\n}"),
                arguments("def ch = flag ? channel.empty() : channel.of(record(id: '1'))", "missing", "Unrecognized property `missing` for type Record {\n    id: String\n}"),
                arguments("def ch = channel.empty().mix(channel.of(record(id: '1')))", "missing", "Unrecognized property `missing` for type Record {\n    id: String\n}"),
                arguments("def ch = channel.empty().mix(channel.value(record(id: '1')))", "missing", "Unrecognized property `missing` for type Record {\n    id: String\n}"),
                arguments("def ch = channel.empty() ; ch = ch.mix(channel.of(record(id: '1')))", "missing", "Unrecognized property `missing` for type Record {\n    id: String\n}"),
                arguments("def ch = channel.empty() ; if( flag ) { ch = channel.of(record(id: '1')) }", "missing", "Unrecognized property `missing` for type Record {\n    id: String\n}"),
                arguments("if( !flag ) { ch = channel.empty() } else { ch = channel.of(record(id: '1')) }", "missing", "Unrecognized property `missing` for type Record {\n    id: String\n}"),
                arguments("if( flag ) { ch = channel.of(record(id: '1')) } else { ch = channel.empty() }", "missing", "Unrecognized property `missing` for type Record {\n    id: String\n}"),
                arguments("def ch = channel.empty() ; ch = channel.of(record(id: '1')) ; ch = channel.of('a')", "id", "Assignment target with type Channel<Record {\n    id: String\n}> cannot be assigned to value with type Channel<String>"),
                arguments("def ch = channel.empty() ; ch = null ; ch = 42", "id", "Assignment target with type Channel<E> cannot be assigned to value with type Integer")
            );
        }

        @ParameterizedTest
        @MethodSource("emptyChannels")
        void inferElementTypeOfEmptyChannel(String source, String field, String error) {
            assertThat(check("""
                def flag = true
                %s
                ch.map { r -> r.%s }
                """.formatted(source, field))).containsExactly(error);
        }

        static Stream<Arguments> emptyValues() {
            return Stream.of(
                arguments("def xs = [] ; if( flag ) { xs = ['a'] } ; xs.each { s -> s.missing }", "Unrecognized property `missing` for type String"),
                arguments("def m = [:] ; if( flag ) { m = [k: 'v'] } ; m.each { k, v -> v.missing }", "Unrecognized property `missing` for type String"),
                arguments("def xs = [] ; xs = ['a'] ; xs = [1]", "Assignment target with type List<String> cannot be assigned to value with type List<Integer>"),
                arguments("def xs = [] ; xs = xs + ['a'] ; xs.each { s -> s.missing }", "Unrecognized property `missing` for type String"),
                arguments("def xs = [] ; xs += ['a'] ; xs.each { s -> s.missing }", "Unrecognized property `missing` for type String"),
                arguments("def e = [] ; def xs = flag ? ['a'] : e ; xs.each { s -> s.missing }", "Unrecognized property `missing` for type String"),
                arguments("def xs = [[]] ; xs = [['a']] ; xs.each { ys -> ys.each { s -> s.missing } }", "Unrecognized property `missing` for type String"),
                arguments("def xs = flag ? [['a']] : [[]] ; xs.each { ys -> ys.each { s -> s.missing } }", "Unrecognized property `missing` for type String"),
                arguments("def t = flag ? tuple(1, []) : tuple(1, ['a']) ; t[1].each { s -> s.missing }", "Unrecognized property `missing` for type String"),
                arguments("def xs = [] ; xs = null ; xs = 'abc'", "Assignment target with type List<E> cannot be assigned to value with type String"),
                arguments("def xs = [] ; xs = [].first() ; xs = 'abc'", "Assignment target with type List<E> cannot be assigned to value with type String"),
                arguments("channel.empty().mix(channel.of(1)).mix(channel.of('a'))", "Function `mix` (with multiple signatures) was called with incorrect number of arguments and/or incorrect argument types")
            );
        }

        @ParameterizedTest
        @MethodSource("emptyValues")
        void inferTypeOfEmptyValue(String source, String error) {
            assertThat(check("""
                def flag = true
                %s
                """.formatted(source))).containsExactly(error);
        }

        @ParameterizedTest
        @ValueSource(strings = {
            "def xs = [] ; xs.each { s -> s.length() }",
            "def m = [:] ; m.each { k, v -> v.length() }",
            "def ys = [] ; def zs = ys[[0]] ; zs.each { z -> z.length() }",
            "channel.topic('x').mix(channel.of(1)).mix(channel.of('a'))",
        })
        void acceptUnknownElementTypeOfEmptyValue(String source) {
            assertThat(getErrors("""
                def flag = true
                %s
                """.formatted(source))).isEmpty();
        }
    }

    @Nested
    class ChannelOperators {

        @Test
        void resolveCombineOfTwoChannels() {
            assertThat(typeOf("""
                left  = channel.of( 42 )
                right = channel.of( 'hello' )
                left.combine(right)
                """)).isEqualTo("Channel<Tuple<Integer, String>>");
        }

        @Test
        void resolveCombineOfTupleChannelsAsFlattenedTuple() {
            assertThat(typeOf("""
                left  = channel.of( tuple(1, 2) )
                right = channel.of( tuple('red', 'blue') )
                left.combine(right)
                """)).isEqualTo("Channel<Tuple<Integer, Integer, String, String>>");
        }

        @Test
        void resolveCombineWithNamedArgumentsAsRecord() {
            assertThat(typeOf("""
                samples = channel.of( record(id: 1, fastq: file('1.fq')) )
                index = channel.value( file('index.fa') )
                samples.combine( single_end: true, index: index )
                """)).isEqualTo("Channel<Record {\n    id: Integer\n    fastq: Path\n    single_end: Boolean\n    index: Path\n}>");
        }

        @Test
        void resolveCombineWithRawTuple() {
            assertThat(typeOf("""
                def pair: Tuple = tuple(1, 'a')
                left  = channel.of( 42 )
                right = channel.of( pair )
                left.combine(right)
                """)).isEqualTo("Channel<Tuple>");
        }

        @Test
        void resolveCombineOfTwoValues() {
            assertThat(typeOf("""
                left  = channel.value( 42 )
                right = channel.value( 'hello' )
                left.combine(right)
                """)).isEqualTo("Value<Tuple<Integer, String>>");
        }

        @Test
        void resolveCombineOfValueWithNamedArgumentsAsRecord() {
            assertThat(typeOf("""
                sample = channel.value( record(id: 1, fastq: file('1.fq')) )
                index = channel.value( file('index.fa') )
                sample.combine( single_end: true, index: index )
                """)).isEqualTo("Value<Record {\n    id: Integer\n    fastq: Path\n    single_end: Boolean\n    index: Path\n}>");
        }

        @Test
        void reportCombineNamedArgumentThatIsChannel() {
            assertThat(check("""
                samples = channel.of( record(id: 1, fastq: file('1.fq')) )
                index = channel.of( file('index.fa') )
                samples.combine( index: index )
                """)).containsExactly("Named argument of a `combine` operation cannot be a channel -- it must be a dataflow value or plain value");
        }

        @ParameterizedTest
        @ValueSource(strings = {
            "left = channel.of( tuple(42, 'hello'), tuple(42, 'goodbye') )\nleft.groupBy()\n",
            "left = channel.of( tuple(42, 2, 'hello'), tuple(42, 2, 'goodbye') )\nleft.groupBy()\n",
        })
        void resolveGroupByOnFirstElementAndBagOfLast(String source) {
            assertThat(typeOf(source)).isEqualTo("Channel<Tuple<Integer, Bag<String>>>");
        }

        @Test
        void resolveJoinOfRecordsAsMergedRecord() {
            assertThat(typeOf("""
                left  = channel.of( record(id: 42, name: 'hello') )
                right = channel.of( record(id: 42, alive: true) )
                left.join(right, by: 'id')
                """)).isEqualTo("Channel<Record {\n    id: Integer\n    name: String\n    alive: Boolean\n}>");
        }

        @Test
        void reportJoinFieldMissingFromBothRecords() {
            assertThat(check("""
                left  = channel.of( record(id: 42, name: 'hello') )
                right = channel.of( record(id: 42, alive: true) )
                left.join(right, by: 'sample_id')
                """)).containsExactly(
                    "Join field `sample_id` is not present in left-hand side",
                    "Join field `sample_id` is not present in right-hand side");
        }

        @Test
        void reportSoftErrorForJoinThatOverwritesLeftHandFields() {
            var errors = getErrors("""
                left  = channel.of( record(id: 42, name: 'hello', alive: false) )
                right = channel.of( record(id: 42, alive: true) )
                left.join(right, by: 'id')
                """);
            assertThat(errors).singleElement().satisfies(e -> {
                assertThat(SeverityAware.isSoftError(e)).isTrue();
                assertThat(e.getOriginalMessage()).isEqualTo("Join fields [alive] are present in both records -- the left-hand values will be overwritten by the right-hand values");
            });
        }
    }

}
