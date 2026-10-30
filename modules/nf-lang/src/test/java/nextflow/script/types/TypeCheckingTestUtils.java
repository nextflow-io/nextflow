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
import org.codehaus.groovy.ast.ClassHelper;
import org.codehaus.groovy.ast.expr.ConstantExpression;
import org.codehaus.groovy.ast.expr.Expression;
import org.codehaus.groovy.ast.stmt.BlockStatement;
import org.codehaus.groovy.ast.stmt.ExpressionStatement;
import org.codehaus.groovy.control.SourceUnit;
import org.codehaus.groovy.control.messages.SyntaxErrorMessage;
import org.codehaus.groovy.syntax.SyntaxException;

import static nextflow.script.types.TypeCheckingUtils.getType;
import static org.assertj.core.api.Assertions.assertThat;

/**
 * Helpers to resolve and type-check a script with static types enabled.
 */
final class TypeCheckingTestUtils {

    private TypeCheckingTestUtils() {}

    /**
     * Parse, resolve and type-check a script with static typing enabled.
     */
    static SourceUnit parse(String contents) {
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
    static List<SyntaxException> getErrors(String contents) {
        var errors = parse(contents).getErrorCollector().getErrors();
        if( errors == null )
            return List.of();
        return errors.stream()
            .filter(e -> e instanceof SyntaxErrorMessage)
            .map(e -> ((SyntaxErrorMessage) e).getCause())
            .sorted(TestUtils.ERROR_COMPARATOR)
            .toList();
    }

    static List<String> check(String contents) {
        return getErrors(contents).stream()
            .map(SyntaxException::getOriginalMessage)
            .toList();
    }

    /**
     * Get the expression of the last statement in the entry workflow
     * of a script that type-checks without errors.
     */
    static Expression parseExpression(String contents) {
        var source = parse(contents);
        assertThat(source.getErrorCollector().getErrors()).isNullOrEmpty();
        var main = ((ScriptNode) source.getAST()).getEntry().main;
        assertThat(main).isInstanceOf(BlockStatement.class);
        var statements = ((BlockStatement) main).getStatements();
        var last = statements.get(statements.size() - 1);
        assertThat(last).isInstanceOf(ExpressionStatement.class);
        return ((ExpressionStatement) last).getExpression();
    }

    /**
     * Get the type name of the last statement in the entry workflow
     * of a script that type-checks without errors.
     */
    static String typeOf(String contents) {
        return Types.getName(getType(parseExpression(contents)));
    }

    static void assertType(Expression exp, Class<?> type) {
        var actual = getType(exp);
        var expected = ClassHelper.makeCached(type);
        assertThat(Types.isEqual(actual, expected))
            .withFailMessage("Expected type of `%s` to be %s but was %s", exp.getText(), Types.getName(expected), Types.getName(actual))
            .isTrue();
    }

}
