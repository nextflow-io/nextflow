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
package nextflow.script.parser;

import java.util.List;

import nextflow.script.ast.ScriptNode;
import nextflow.script.control.ScriptParser;
import nextflow.script.control.ScriptToGroovyHelper;
import nextflow.test.TestUtils;
import org.codehaus.groovy.ast.expr.ConstantExpression;
import org.codehaus.groovy.ast.expr.ConstructorCallExpression;
import org.codehaus.groovy.ast.expr.ListExpression;
import org.codehaus.groovy.ast.expr.TupleExpression;
import org.codehaus.groovy.control.SourceUnit;
import org.codehaus.groovy.syntax.SyntaxException;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;
import org.junit.jupiter.params.provider.ValueSource;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.tuple;

/**
 * Specifies the parsing and analysis of agent definitions.
 *
 * @see ScriptAstBuilder
 */
class AgentParserTest {

    private static List<SyntaxException> check(String contents) {
        return TestUtils.check(new ScriptParser(), contents);
    }

    private static SourceUnit parse(String contents) {
        var parser = new ScriptParser();
        var source = parser.parse("main.nf", contents);
        parser.analyze();
        assertThat(TestUtils.hasSyntaxErrors(source)).isFalse();
        return source;
    }

    private static ScriptNode ast(SourceUnit source) {
        return (ScriptNode) source.getAST();
    }

    private static List<String> messages(List<SyntaxException> errors) {
        return errors.stream().map(SyntaxException::getOriginalMessage).toList();
    }

    // The prompt closure's free-variable refs (params.*, task.ext.*) are captured via the same
    // collector that populates process-body BodyDef.valRefs, so AgentToGroovyVisitor can fold
    // them into the synthetic PromptDef/BodyDef cache key.

    @Test
    void capturePromptGlobalsAsValRefsExcludingDeclaredInputs() {
        var source = parse("""
            nextflow.enable.types = true

            agent eval_agent {
                model 'openai/gpt-5-mini'
                tools()

                input:
                question: String

                output:
                stdout()

                prompt:
                \"""
                ${question} threshold=${params.threshold} args=${task.ext.args}
                \"""
            }
            """);
        var node = ast(source).getAgents().get(0);

        var refs = (ListExpression) new ScriptToGroovyHelper(source).getVariableRefs(node.prompt);
        var names = refs.getExpressions().stream()
            .map(ref -> ((TupleExpression) ((ConstructorCallExpression) ref).getArguments()).getExpression(0))
            .map(name -> ((ConstantExpression) name).getValue())
            .toList();

        // the declared input `question` is not captured
        assertThat(names).containsExactlyInAnyOrder("params.threshold", "task.ext.args");
    }

    @Test
    void captureNoPromptValRefsWhenPromptReferencesOnlyDeclaredInputs() {
        var source = parse("""
            nextflow.enable.types = true

            agent eval_agent {
                model 'openai/gpt-5-mini'
                tools()

                input:
                question: String

                output:
                stdout()

                prompt:
                \"""
                Question: ${question}
                \"""
            }
            """);
        var node = ast(source).getAgents().get(0);

        var refs = (ListExpression) new ScriptToGroovyHelper(source).getVariableRefs(node.prompt);

        assertThat(refs.getExpressions()).isEmpty();
    }

    @Test
    void parseMinimalAgentDefinition() {
        var script = ast(parse("""
            nextflow.enable.types = true

            record Question {
                text: String
            }

            record Answer {
                plan: String
            }

            agent eval_agent {
                model 'openai/gpt-5-mini'
                instruction 'You are helpful.'
                tools()
                maxIterations 20

                input:
                question: Question

                output:
                plan: Answer

                prompt:
                \"""
                Question: ${question.text}
                \"""
            }
            """));

        assertThat(script.getAgents()).hasSize(1);
        var node = script.getAgents().get(0);
        assertThat(node.getName()).isEqualTo("eval_agent");
        assertThat(node.inputs).hasSize(1);
        assertThat(node.inputs[0].getName()).isEqualTo("question");
    }

    @Test
    void reportAgentWithoutPromptSection() {
        var errors = check("""
            nextflow.enable.types = true

            record Question {
                text: String
            }

            agent broken {
                model 'openai/gpt-5-mini'
                instruction 'x'
                tools()

                input:
                q: Question
            }
            """);

        assertThat(errors)
            .extracting(SyntaxException::getStartLine, SyntaxException::getOriginalMessage)
            .containsExactly(tuple(14, "Unexpected input: '}'"));
    }

    @Test
    void resolveAgentReferenceFromWorkflow() {
        var errors = check("""
            nextflow.enable.types = true

            record Question {
                text: String
            }

            record Answer {
                text: String
            }

            agent eval_agent {
                model 'openai/gpt-5-mini'
                instruction 'x'
                tools()

                input:
                q: Question

                output:
                r: Answer

                prompt:
                \"""
                ${q.text}
                \"""
            }

            workflow {
                answer = eval_agent(record(text: 'hi'))
                answer.view()
            }
            """);

        assertThat(errors).isEmpty();
    }

    @Test
    void reportAgentUsedAsVariable() {
        var errors = check("""
            nextflow.enable.types = true

            agent eval_agent {
                model 'openai/gpt-5-mini'
                instruction 'x'

                input:
                q: String

                output:
                stdout()

                prompt:
                \"""
                ${q}
                \"""
            }

            workflow {
                channel.of('hi') | eval_agent
            }
            """);

        assertThat(errors)
            .extracting(SyntaxException::getStartLine, SyntaxException::getStartColumn, SyntaxException::getOriginalMessage)
            .containsExactly(tuple(20, 24, "Agent `eval_agent` cannot be used as a variable"));
    }

    @Test
    void resolveDirectivesAndPromptVariablesInAgentBody() {
        var errors = check("""
            nextflow.enable.types = true

            record Question {
                text: String
            }

            record Answer {
                plan: String
            }

            agent eval_agent {
                model 'openai/gpt-5-mini'
                instruction 'You are helpful.'
                tools()
                maxIterations 20

                input:
                question: Question

                output:
                plan: Answer

                prompt:
                \"""
                Question: ${question.text}
                \"""
            }
            """);

        assertThat(errors).isEmpty();
    }

    @Test
    void resolveRepeatableLabelDirectiveInAgentBody() {
        var errors = check("""
            nextflow.enable.types = true

            agent eval_agent {
                label 'reasoning'
                label 'fast'
                model 'openai/gpt-5-mini'

                input:
                question: String

                output:
                stdout()

                prompt:
                \"""
                Question: ${question}
                \"""
            }
            """);

        assertThat(errors).isEmpty();
    }

    @Test
    void acceptRecordTypedAgentIO() {
        var errors = check("""
            nextflow.enable.types = true

            record Question {
                text: String
            }

            record Answer {
                answer: String
            }

            agent eval_agent {
                model 'openai/gpt-5-mini'
                instruction 'You are helpful.'
                tools()

                input:
                q: Question

                output:
                a: Answer

                prompt:
                \"""
                ${q.text}
                \"""
            }
            """);

        assertThat(errors).isEmpty();
    }

    @Test
    void acceptValAgentIO() {
        var errors = check("""
            nextflow.enable.types = true

            agent eval_agent {
                model 'openai/gpt-5-mini'
                instruction 'You are helpful.'
                tools()

                input:
                question: String

                output:
                stdout()

                prompt:
                \"""
                ${question}
                \"""
            }
            """);

        assertThat(errors).isEmpty();
    }

    @Test
    void acceptDestructuredAgentInputs() {
        var errors = check("""
            nextflow.enable.types = true

            agent qa {
                input:
                record(id: String, reads: Path)
                tuple(a: Integer, b: Path)

                output:
                record(id: id, a: a, summary: stdout())

                prompt:
                "Summarize ${id} ${reads} ${a} ${b}"
            }
            """);

        assertThat(errors).isEmpty();
    }

    @Test
    void acceptComposedAgentOutput() {
        var errors = check("""
            nextflow.enable.types = true

            agent qa {
                input:
                q: String

                output:
                file('report.md')

                prompt:
                "go"
            }
            """);

        assertThat(errors).isEmpty();
    }

    @Test
    void rejectMultipleAgentOutputs() {
        var errors = check("""
            nextflow.enable.types = true

            agent qa {
                input:
                q: String

                output:
                count: Integer
                report: Path = file('report.md')

                prompt:
                "go"
            }
            """);

        assertThat(messages(errors)).contains("Agent should have only one output -- combine outputs into a record");
    }

    @Test
    void allowNamedOutputWhenItIsTheOnlyAgentOutput() {
        var errors = check("""
            nextflow.enable.types = true

            agent qa {
                input:
                q: String

                output:
                summary: String = stdout()

                prompt:
                "go"
            }
            """);

        assertThat(errors).isEmpty();
    }

    private static List<String> agentOutputErrors(String declaration) {
        var errors = check("""
            nextflow.enable.types = true

            record Answer {
                text: String
            }

            record BadAnswer {
                text: String
                meta: Map
            }

            agent qa {
                input:
                q: String

                output:
                %s

                prompt:
                "go"
            }
            """.formatted(declaration));
        return messages(errors).stream()
            .filter(message -> message.startsWith("Agent output"))
            .toList();
    }

    @ParameterizedTest
    @ValueSource(strings = {
        "a: Answer",
        "n: Integer",
        "x: Float",
        "ok: Boolean",
        "f: Path",
        "s: String",
        "items: List<String>",
        "answers: List<Answer>",
        "nested: List<List<Path>>",
    })
    void acceptSupportedModelAnsweredOutputTypes(String declaration) {
        assertThat(agentOutputErrors(declaration)).isEmpty();
    }

    @ParameterizedTest
    @CsvSource(delimiter = '|', textBlock = """
        items: List         | Agent output `items` has unsupported type List -- supported types are Boolean, Float, Integer, List<E>, Path, String, or a record type
        items: List<Map>    | Agent output `items` has unsupported type List<Map> -- supported types are Boolean, Float, Integer, List<E>, Path, String, or a record type
        d: Double           | Agent output `d` has unsupported type Double -- supported types are Boolean, Float, Integer, List<E>, Path, String, or a record type
        m: Map              | Agent output `m` has unsupported type Map -- supported types are Boolean, Float, Integer, List<E>, Path, String, or a record type
        b: BadAnswer        | Agent output `b` has unsupported field `meta` with type Map -- supported types are Boolean, Float, Integer, List<E>, Path, String, or a record type
        bs: List<BadAnswer> | Agent output `bs` has unsupported field `meta` with type Map -- supported types are Boolean, Float, Integer, List<E>, Path, String, or a record type
        q                   | Agent output `q` should declare a type -- typed outputs are answered by the model
        """)
    void rejectUnsupportedModelAnsweredOutputTypes(String declaration, String expected) {
        assertThat(agentOutputErrors(declaration)).containsExactly(expected);
    }

    @Test
    void resolveFileFilesAndStdoutInAgentOutput() {
        var errors = check("""
            nextflow.enable.types = true

            agent qa {
                input:
                q: String

                output:
                record(
                    summary: stdout(),
                    report: file('report.md'),
                    notes: files('*.txt')
                )

                prompt:
                "go"
            }
            """);

        assertThat(errors).isEmpty();
    }

    @Test
    void doNotResolveProcessOnlyEvalInAgentOutput() {
        var errors = check("""
            nextflow.enable.types = true

            agent qa {
                input:
                q: String

                output:
                eval('date')

                prompt:
                "go"
            }
            """);

        assertThat(messages(errors)).anySatisfy(message -> assertThat(message).contains("eval"));
    }

}
