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

import java.nio.file.Path;
import java.util.List;
import java.util.Map;
import java.util.stream.Stream;

import nextflow.script.ast.ScriptNode;
import nextflow.script.control.ScriptParser;
import nextflow.script.namespaces.ChannelNamespace;
import org.codehaus.groovy.ast.ClassHelper;
import org.codehaus.groovy.ast.ClassNode;
import org.codehaus.groovy.ast.GenericsType;
import org.codehaus.groovy.ast.GenericsType.GenericsTypeName;
import org.codehaus.groovy.ast.expr.Expression;
import org.codehaus.groovy.ast.stmt.BlockStatement;
import org.codehaus.groovy.ast.stmt.ExpressionStatement;
import org.codehaus.groovy.ast.tools.GenericsUtils;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.Arguments;
import org.junit.jupiter.params.provider.MethodSource;

import static nextflow.script.types.TypeCheckingTestUtils.assertType;
import static nextflow.script.types.TypeCheckingUtils.*;
import static org.assertj.core.api.Assertions.assertThat;
import static org.junit.jupiter.params.provider.Arguments.arguments;

/**
 * Specifies the type inference and generics helpers used by the type checker.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class TypeCheckingUtilsTest {

    private static final ClassNode CHANNEL_TYPE = ClassHelper.makeCached(Channel.class);

    private static final GenericsTypeName E = new GenericsTypeName("E");

    /**
     * Parse and analyze a script and return the expression of the
     * last statement in the entry workflow.
     */
    private static Expression parseExpression(String contents) {
        var parser = new ScriptParser();
        var source = parser.parse("main.nf", contents);
        parser.analyze();
        assertThat(source.getErrorCollector().getErrors()).isNullOrEmpty();
        var main = ((ScriptNode) source.getAST()).getEntry().main;
        assertThat(main).isInstanceOf(BlockStatement.class);
        var statements = ((BlockStatement) main).getStatements();
        var last = statements.get(statements.size() - 1);
        assertThat(last).isInstanceOf(ExpressionStatement.class);
        return ((ExpressionStatement) last).getExpression();
    }

    static Stream<Arguments> literals() {
        return Stream.of(
            arguments("true", Boolean.class),
            arguments("42", Integer.class),
            arguments("3.14", Float.class),
            arguments("'hello'", String.class),
            arguments("\"hello ${\"world!\"}\"", String.class),
            arguments("[1, 2, 3]", List.class),
            arguments("[a: 1, b: 2, c: 3]", Map.class)
        );
    }

    @ParameterizedTest
    @MethodSource("literals")
    void resolveTypeOfLiteral(String expr, Class<?> type) {
        assertType(parseExpression(expr), type);
    }

    @Test
    void resolveReturnTypeOfFunctionCall() {
        assertType(parseExpression("env('HELLO')"), String.class);
    }

    static Stream<Arguments> propertyAccesses() {
        return Stream.of(
            arguments("workflow.projectDir", Path.class),
            arguments("workflow.projectDir.name", String.class)
        );
    }

    @ParameterizedTest
    @MethodSource("propertyAccesses")
    void resolveTypeOfPropertyAccess(String expr, Class<?> type) {
        assertType(parseExpression(expr), type);
    }

    @Test
    void resolveReturnTypeOfMethodCall() {
        assertType(parseExpression("workflow.projectDir.exists()"), Boolean.class);
    }

    @Test
    void resolveOutputTypeOfProcessInvocation() {
        var exp = parseExpression("""
            nextflow.enable.types = true

            process hello {
                input:
                target: String

                output:
                x = 'hello'
                y = target

                script:
                \"""
                \"""
            }

            workflow {
                hello('world')
            }
            """);
        var cn = getType(exp);
        assertThat(cn.getField("x")).isNotNull();
        assertThat(cn.getField("y")).isNotNull();
    }

    @Test
    void resolveOutputTypeOfWorkflowInvocation() {
        var exp = parseExpression("""
            nextflow.enable.types = true

            workflow hello {
                take:
                target: String

                emit:
                x = 'hello'
                y = target
            }

            workflow {
                hello('world')
            }
            """);
        var cn = getType(exp);
        assertThat(cn.getField("x")).isNotNull();
        assertThat(cn.getField("y")).isNotNull();
    }

    @Test
    void extractPlaceholderMappingFromParameterizedType() {
        assertThat(GenericsUtils.extractPlaceholders(CHANNEL_TYPE)).isEmpty();

        var spec = GenericsUtils.extractPlaceholders(makeType(CHANNEL_TYPE, ClassHelper.Integer_TYPE));
        assertThat(spec.get(E)).extracting(GenericsType::getType).isEqualTo(ClassHelper.Integer_TYPE);
    }

    @Test
    void applyPlaceholderMappingToParameterizedType() {
        var declaringType = makeType(CHANNEL_TYPE, ClassHelper.Integer_TYPE);

        var cn = resolveGenericType(CHANNEL_TYPE, declaringType);
        assertThat(GenericsUtils.extractPlaceholders(cn).get(E))
            .extracting(GenericsType::getType)
            .isEqualTo(ClassHelper.Integer_TYPE);

        var mn = ClassHelper.makeCached(ChannelNamespace.class).getDeclaredMethods("of").get(0);
        var returnType = resolveGenericType(mn.getReturnType(), declaringType);
        assertThat(GenericsUtils.extractPlaceholders(returnType).get(E))
            .extracting(GenericsType::getType)
            .isEqualTo(ClassHelper.Integer_TYPE);
    }

}
