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

import java.util.Collections;
import java.util.List;
import java.util.stream.Stream;

import nextflow.script.ast.ASTNodeMarker;
import nextflow.script.ast.AssignmentExpression;
import nextflow.script.ast.ProcessNodeV2;
import nextflow.script.ast.ScriptNode;
import nextflow.test.TestUtils;
import org.codehaus.groovy.ast.MethodNode;
import org.codehaus.groovy.ast.expr.ClosureExpression;
import org.codehaus.groovy.ast.expr.ConstantExpression;
import org.codehaus.groovy.ast.expr.Expression;
import org.codehaus.groovy.ast.expr.MethodCallExpression;
import org.codehaus.groovy.ast.expr.PropertyExpression;
import org.codehaus.groovy.ast.expr.TupleExpression;
import org.codehaus.groovy.ast.expr.VariableExpression;
import org.codehaus.groovy.ast.stmt.BlockStatement;
import org.codehaus.groovy.ast.stmt.ExpressionStatement;
import org.codehaus.groovy.ast.stmt.Statement;
import org.codehaus.groovy.control.SourceUnit;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.Arguments;
import org.junit.jupiter.params.provider.CsvSource;
import org.junit.jupiter.params.provider.MethodSource;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.tuple;
import static org.junit.jupiter.params.provider.Arguments.arguments;

/**
 * Specifies the lowering of agents by {@link AgentToGroovyVisitor}.
 *
 * An agent's typed I/O lowers through the same compiler units a process's does, so the implicit
 * stagers and the output unstagers must come out identical for identical declarations.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AgentToGroovyTest {

    static Stream<Arguments> pathTypedInputs() {
        return Stream.of(
            arguments("contigs: Path", List.of("contigs")),
            arguments("contigs: Path?", List.of("contigs")),
            arguments("reads: Set<Path>", List.of("reads")),
            arguments("reads: List<Path>", List.of("reads")),
            arguments("n: Integer", List.of()),
            arguments("s: String", List.of()),
            arguments("contigs: Path\nn: Integer", List.of("contigs")),
            arguments("record(id: String, reads: Path)", List.of("reads")),
            arguments("tuple(id: String, reads: Path)", List.of("reads"))
        );
    }

    @ParameterizedTest
    @MethodSource("pathTypedInputs")
    void generateImplicitStagerForEachPathTypedInput(String declaration, List<String> expected) {
        var blocks = lowerAgent("""
            agent qa {
                input:
                %s
                output:
                stdout()
                prompt: "go"
            }
            """.formatted(declaration));
        assertThat(stagerTargets(blocks.stagers())).isEqualTo(expected);
    }

    @Test
    void lowerDestructuredInputIntoTupleInput() {
        var blocks = lowerAgent("""
            agent qa {
                input:
                record(id: String, reads: Path)
                output:
                record(id: id, summary: stdout())
                prompt: "go"
            }
            """);
        var input = firstCall(blocks.inputs());
        var output = firstCall(blocks.outputs());

        assertThat(input.getMethodAsString()).isEqualTo("_input_");
        assertThat(args(input)).hasSize(2);
        assertThat(args(input).get(1).getType().getName()).isEqualTo("nextflow.script.types.Record");

        // a bare expression becomes the implicit output, with stdout() decoded to the answer
        assertThat(args(output).get(0).getText()).isEqualTo("$out");
        assertThat(closureBody((ClosureExpression) args(output).get(2)).getText())
            .contains("nextflow.agent.AgentOutputPlan.answer(this.stdout())");
    }

    @ParameterizedTest
    @CsvSource(delimiter = '|', textBlock = """
        stdout()          | java.lang.String
        file('report.md') | java.nio.file.Path
        """)
    void useInferredTypeForUnnamedOutput(String declaration, String expected) {
        var blocks = lowerAgent("""
            agent qa {
                input:
                q: String
                output:
                %s
                prompt: "go"
            }
            """.formatted(declaration));
        var output = firstCall(blocks.outputs());

        assertThat(args(output).get(1).getType().getName()).isEqualTo(expected);
    }

    @Test
    void preserveParameterizedOutputTypeInHiddenClass() {
        var blocks = lowerAgent("""
            agent qa {
                input:
                q: String
                output:
                scores: List<Integer>
                prompt: "go"
            }
            """);
        var args = args(firstCall(blocks.outputs()));

        assertThat(args).hasSize(3);
        assertThat(args.get(0).getText()).isEqualTo("scores");
        assertThat(args.get(2).getText()).isEqualTo("scores");
        assertThat(args.get(1).getType().getField("scores").getType().toString(false))
            .isEqualTo("java.util.List<java.lang.Integer>");
    }

    @Test
    void generateSameStagersAsEquivalentProcess() {
        var source = analyze("""
            nextflow.enable.types = true

            record Sample {
                id: String
                seq: Path
                index: Path
            }

            process CHECK {
                input:
                sample: Sample
                n: Integer
                output:
                'x'
                script:
                'true'
            }

            agent qa {
                input:
                sample: Sample
                n: Integer
                output:
                stdout()
                prompt: "go"
            }
            """);
        var ast = (ScriptNode) source.getAST();

        // the record recursion runs for both
        var agent = new AgentToGroovyVisitor(source).transform(ast.getAgents().get(0));
        var process = new ProcessToGroovyVisitorV2(source).transform((ProcessNodeV2) ast.getProcesses().get(0));

        assertThat(stagerTargets(blockAt(agent, 1))).containsExactly("sample.seq", "sample.index");
        assertThat(stagerTargets(blockAt(process, 1))).isEqualTo(stagerTargets(blockAt(agent, 1)));
    }

    @Test
    void carryNullableMarkerIntoInputDeclaration() {
        var blocks = lowerAgent("""
            agent qa {
                input:
                a: Path
                b: Path?
                output:
                stdout()
                prompt: "go"
            }
            """);

        // the third argument of _input_ is the optional flag
        assertThat(calls(blocks.inputs()))
            .extracting(
                call -> ((ConstantExpression) args(call).get(0)).getValue(),
                call -> args(call).get(1).getType().getName(),
                call -> ((ConstantExpression) args(call).get(2)).getValue())
            .containsExactly(
                tuple("a", "java.nio.file.Path", false),
                tuple("b", "java.nio.file.Path", true));
    }

    @Test
    void lowerFileOutputIntoUnstagerPlusValueCarryingOutput() {
        var blocks = lowerAgent("""
            agent qa {
                input:
                q: String
                output:
                record(
                    report: file('report.md'),
                    notes: files('*.txt')
                )
                prompt: "go"
            }
            """);

        // one unstager per file()/files() call, keyed per agent
        assertThat(calls(blocks.unstagers()))
            .extracting(
                MethodCallExpression::getMethodAsString,
                call -> ((ConstantExpression) args(call).get(0)).getValue())
            .containsExactly(
                tuple("_unstage_files", "$path0"),
                tuple("_unstage_files", "$path1"));

        // the output carries its closure
        assertThat(calls(blocks.outputs()))
            .extracting(call -> args(call).size())
            .containsExactly(3);
    }

    @Test
    void doNotLowerEnvOrEvalInOutput() {
        var blocks = lowerAgent("""
            agent qa {
                input:
                q: String
                output:
                env('HOME')
                prompt: "go"
            }
            """);

        // an agent has no task environment to read back
        assertThat(blocks.unstagers().getStatements()).isEmpty();
    }

    @Test
    void placeStagersAndUnstagersBeforePrompt() {
        // the prompt must stay last: its value is the closure return value
        var blocks = lowerAgent("""
            agent qa {
                input:
                contigs: Path
                output:
                file('report.md')
                prompt: "go"
            }
            """);

        assertThat(stagerTargets(blocks.stagers())).containsExactly("contigs");
        assertThat(calls(blocks.unstagers()))
            .extracting(MethodCallExpression::getMethodAsString)
            .containsExactly("_unstage_files");
    }

    @Test
    void lowerEveryFileCallFormProcessOutputAccepts() {
        // the one-arg and two-arg opts forms, singular and plural. DSL resolution matches by
        // name only, so this covers the lowering rather than the declared signature
        var blocks = lowerAgent("""
            nextflow.enable.types = true

            agent qa {
                input:
                x: String
                output:
                record(
                    one: file('report.md'),
                    typed: file(type: 'file', 'report.md'),
                    many: files('*.tsv'),
                    opts: files(hidden: true, '*.log')
                )
                prompt: "go"
            }
            """);

        // all four resolve, and each lowers to its own unstager
        assertThat(calls(blocks.unstagers()))
            .extracting(MethodCallExpression::getMethodAsString)
            .isEqualTo(Collections.nCopies(4, "_unstage_files"));
    }

    @Test
    void resolveFileInOutputToSameScopeAsProcessOutput() {
        var source = analyze("""
            nextflow.enable.types = true

            process CHECK {
                input:
                x: String
                output:
                file('report.md')
                script:
                "true"
            }

            agent qa {
                input:
                x: String
                output:
                file('report.md')
                prompt: "go"
            }
            """);
        var ast = (ScriptNode) source.getAST();

        var agentTarget = fileCallTarget(ast.getAgents().get(0).outputs);
        var processTarget = fileCallTarget(((ProcessNodeV2) ast.getProcesses().get(0)).outputs);

        // the shared declaration, not the driver-side global `ScriptDsl.file`, which also matches by name
        assertThat(agentTarget).isEqualTo("nextflow.script.dsl.FileOutputDsl");
        // provably one declaration, not two that may drift
        assertThat(processTarget).isEqualTo(agentTarget);
    }

    private static SourceUnit analyze(String contents) {
        var parser = new ScriptParser();
        var source = parser.parse("main.nf", contents);
        parser.analyze();
        assertThat(TestUtils.getErrors(source)).isEmpty();
        return source;
    }

    /**
     * The blocks of the generated agent closure body, in order:
     * directives, stagers, unstagers, inputs, outputs, prompt.
     */
    private record LoweredAgent(
        BlockStatement stagers,
        BlockStatement unstagers,
        BlockStatement inputs,
        BlockStatement outputs
    ) {}

    private static LoweredAgent lowerAgent(String contents) {
        var source = analyze(contents);
        var agent = ((ScriptNode) source.getAST()).getAgents().get(0);
        var stmt = new AgentToGroovyVisitor(source).transform(agent);
        return new LoweredAgent(
            blockAt(stmt, 1),
            blockAt(stmt, 2),
            blockAt(stmt, 3),
            blockAt(stmt, 4));
    }

    /**
     * The i-th statement of the generated `agent('name') { ... }` closure body.
     */
    private static BlockStatement blockAt(Statement stmt, int index) {
        var call = (MethodCallExpression) ((ExpressionStatement) stmt).getExpression();
        var closure = (ClosureExpression) args(call).get(1);
        return (BlockStatement) ((BlockStatement) closure.getCode()).getStatements().get(index);
    }

    private static List<MethodCallExpression> calls(BlockStatement block) {
        return block.getStatements().stream()
            .map(st -> (MethodCallExpression) ((ExpressionStatement) st).getExpression())
            .toList();
    }

    private static MethodCallExpression firstCall(BlockStatement block) {
        return calls(block).get(0);
    }

    private static List<Expression> args(MethodCallExpression call) {
        return ((TupleExpression) call.getArguments()).getExpressions();
    }

    /**
     * The `stageAs({ <target> })` targets, rendered as source-like text.
     */
    private static List<String> stagerTargets(BlockStatement block) {
        return calls(block).stream()
            .map(call -> {
                assertThat(call.getMethodAsString()).isEqualTo("stageAs");
                return render(closureBody((ClosureExpression) args(call).get(0)));
            })
            .toList();
    }

    /**
     * The single expression of a one-statement closure, wrapped in a block or not.
     */
    private static Expression closureBody(ClosureExpression closure) {
        var code = closure.getCode();
        var st = code instanceof BlockStatement block ? block.getStatements().get(0) : code;
        return ((ExpressionStatement) st).getExpression();
    }

    private static String render(Expression expr) {
        if( expr instanceof VariableExpression ve )
            return ve.getName();
        if( expr instanceof PropertyExpression pe )
            return render(pe.getObjectExpression()) + "." + pe.getPropertyAsString();
        return expr.getText();
    }

    /**
     * The declaring class the `file(...)` call in the first output declaration resolved to.
     * `file` is overloaded, so the scope visitor records METHOD_OVERLOADS rather than a single
     * METHOD_TARGET; both are read so the assertion does not depend on the arity used.
     */
    @SuppressWarnings("unchecked")
    private static String fileCallTarget(Statement outputs) {
        var st = ((BlockStatement) outputs).getStatements().get(0);
        var expr = ((ExpressionStatement) st).getExpression();
        var call = (MethodCallExpression) (expr instanceof AssignmentExpression ae ? ae.getRightExpression() : expr);
        var overloads = (List<MethodNode>) call.getNodeMetaData(ASTNodeMarker.METHOD_OVERLOADS);
        var mn = overloads != null && !overloads.isEmpty()
            ? overloads.get(0)
            : (MethodNode) call.getNodeMetaData(ASTNodeMarker.METHOD_TARGET);
        if( mn == null || mn.getDeclaringClass() == null )
            return null;
        return mn.getDeclaringClass().getName();
    }

}
