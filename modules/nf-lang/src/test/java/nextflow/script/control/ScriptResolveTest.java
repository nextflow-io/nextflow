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

import java.nio.file.Path;
import java.util.List;

import nextflow.test.TestUtils;
import org.codehaus.groovy.control.SourceUnit;
import org.codehaus.groovy.syntax.SyntaxException;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;

import static nextflow.test.TestUtils.assertErrors;
import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.tuple;

/**
 * Specifies how {@link ScriptResolveVisitor} resolves names in a script.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class ScriptResolveTest {

    private static List<SyntaxException> check(String contents) {
        return TestUtils.check(new ScriptParser(), contents);
    }

    private static List<SyntaxException> check(Path projectDir, Path... files) {
        return TestUtils.check(new ScriptParser(projectDir), List.of(files));
    }

    private static SourceUnit analyze(String contents) {
        var parser = new ScriptParser();
        var source = parser.parse("main.nf", contents);
        parser.analyze();
        return source;
    }

    @Nested
    class ConflictingDeclarations {

        @Test
        void functionAndWorkflow() {
            assertErrors(check("""
                def hello() {
                }

                workflow hello {
                }
                """),
                tuple(4, 1, "`hello` is already declared"));
        }

        @Test
        void functionParamAndLocalVariable() {
            assertErrors(check("""
                def hello(x) {
                    def x = 'world'
                }
                """),
                tuple(2, 9, "`x` is already declared"));
        }

        @Test
        void workflowTakeAndLocalVariable() {
            assertErrors(check("""
                workflow hello {
                    take:
                    x

                    main:
                    def x = 'world'
                }
                """),
                tuple(6, 9, "`x` is already declared"));
        }
    }

    @Nested
    class ProcessSections {

        /**
         * {@code file()} reaches a process output through OutputDslV2 -> FileOutputDsl.
         * Resolution used to follow {@code getInterfaces()[0]} only, so a later
         * {@code extends A, FileOutputDsl} would have silently broken {@code file()}
         * in every process output block.
         */
        @Test
        void resolveDslFunctionInheritedFromSupertypeThatIsNotFirst() {
            assertThat(check("""
                nextflow.enable.types = true

                process CHECK {
                    input:
                    x: String
                    output:
                    file('report.md')
                    script:
                    "true"
                }
                """)).isEmpty();
        }

        /**
         * {@code file}/{@code files} come from the shared FileOutputDsl rather than
         * ProcessDsl.OutputDslV2. This pins that the process-only {@code env},
         * {@code eval} and {@code stdout} survived the extraction. Which {@code file}
         * resolves is pinned in AgentToGroovyTest, because a global
         * {@code ScriptDsl.file} matches by name and would mask a broken scope here.
         */
        @Test
        void resolveEveryTypedProcessOutputDirective() {
            assertThat(check("""
                nextflow.enable.types = true

                process FOO {
                    input:
                    x: String

                    output:
                    report: Path = file('report.md')
                    notes: Set<Path> = files(hidden: true, '*.txt')
                    home: String = env('HOME')
                    lines: String = eval('wc -l report.md')
                    log: String = stdout()

                    script:
                    "true"
                }
                """)).isEmpty();
        }

        @Test
        void reportUnrecognizedDirectiveAndQualifiers() {
            assertErrors(check("""
                process hello {
                    cpu 8

                    input:
                    stdout

                    output:
                    stdin

                    script:
                    ""
                }
                """),
                tuple(2, 5, "Unrecognized process directive `cpu`"),
                tuple(5, 5, "Unrecognized process input qualifier `stdout`"),
                tuple(8, 5, "Unrecognized process output qualifier `stdin`"));
        }

        @Test
        void allowStdinAndStdoutWithoutParens() {
            assertThat(check("""
                process hello {
                    input:
                    tuple val(id), stdin

                    output:
                    tuple val(id), stdout

                    script:
                    ""
                }
                """)).isEmpty();
        }
    }

    @Nested
    class UndefinedVariables {

        @Test
        void inStatement() {
            assertErrors(check("""
                println(x)
                """),
                tuple(1, 9, "`x` is not defined"));
        }

        @Test
        void assignedButNotDeclared() {
            assertErrors(check("""
                def hello() {
                    x = 'hello'
                }
                """),
                tuple(2, 5, "`x` was assigned but not declared"));
        }

        @Test
        void assignedInClosure() {
            assertErrors(check("""
                workflow {
                    [1, 2, 3].each { i ->
                        x = 'hello'
                    }
                }
                """),
                tuple(3, 9, "Variables in a closure should be declared with `def`"));
        }

        @Test
        void inWorkflowEmit() {
            assertErrors(check("""
                workflow hello {
                    emit:
                    x
                }
                """),
                tuple(3, 5, "`x` is not defined"));
        }

        @Test
        void inWorkflowEmitAssignment() {
            assertErrors(check("""
                workflow hello {
                    emit:
                    x = y
                }
                """),
                tuple(3, 9, "`y` is not defined"));
        }

        @Test
        void inWorkflowEmitExpression() {
            assertErrors(check("""
                workflow hello {
                    emit:
                    x + y
                }
                """),
                tuple(3, 5, "`x` is not defined"),
                tuple(3, 9, "`y` is not defined"));
        }

        @Test
        void inTypedWorkflowEmit() {
            assertErrors(check("""
                nextflow.enable.types = true

                workflow hello {
                    emit:
                    x: Channel<String>
                }
                """),
                tuple(5, 5, "`x` is not defined"));
        }
    }

    @Nested
    class Types {

        @Test
        void undefinedTypeInCast() {
            assertErrors(check("""
                'hello' as FooBar
                """),
                tuple(1, 1, "`FooBar` is not defined"));
        }

        @Test
        void undefinedTypeArgument() {
            assertErrors(check("""
                [] as List<FooBar>
                """),
                tuple(1, 1, "`FooBar` is not defined"));
        }

        @Test
        void undefinedTypeArguments() {
            assertErrors(check("""
                [:] as Map<Foo,Bar>
                """),
                tuple(1, 1, "`Foo` is not defined"),
                tuple(1, 1, "`Bar` is not defined"));
        }

        @Test
        void undefinedTypeInConstructorCall() {
            assertErrors(check("""
                x = new FooBar()
                """),
                tuple(1, 5, "`FooBar` is not defined"));
        }

        // requires groovy-json on the test runtime classpath
        @ParameterizedTest
        @ValueSource(strings = {
            "new groovy.json.JsonSlurper()",
            "new groovy.json.JsonGenerator.Options()"
        })
        void recognizeFullyQualifiedClassName(String contents) {
            assertThat(check(contents + "\n")).isEmpty();
        }
    }

    @Nested
    class Warnings {

        @Test
        void warnAboutProcessWhenSection() {
            var source = analyze("""
                nextflow.enable.types = true

                process hello {
                    when:
                    task.ext.when

                    exec:
                    println 'hello!'
                }
                """);
            assertThat(TestUtils.getErrors(source)).isEmpty();
            assertThat(TestUtils.getWarnings(source)).singleElement().asString()
                .contains("The `when` section is deprecated");
        }

        @Test
        void warnWhenPublishStatementHasTypeAnnotation() {
            var source = analyze("""
                nextflow.enable.types = true

                workflow {
                    main:
                    ch = channel.of(1)

                    publish:
                    a: Channel<Integer> = ch
                }

                output {
                    a: Channel<Integer> {
                        path 'a'
                    }
                }
                """);
            assertThat(TestUtils.getErrors(source)).isEmpty();
            assertThat(TestUtils.getWarnings(source)).singleElement().asString()
                .contains("Workflow output type should be declared in the output block");
        }

        @Test
        void doNotWarnWhenWorkflowEmitsChannelByName() {
            var source = analyze("""
                workflow example {
                    main:
                    output = channel.value('hi')

                    emit:
                    output
                }
                """);
            assertThat(TestUtils.getErrors(source)).isEmpty();
            assertThat(TestUtils.getWarnings(source)).isEmpty();
        }

        @Test
        void doNotWarnWhenTypedOutputNamesVariableFromBody() {
            var source = analyze("""
                nextflow.enable.types = true

                process foo {
                    output:
                    message: String

                    exec:
                    message = 'hi'
                }

                workflow example {
                    main:
                    result = foo()

                    emit:
                    result: Value<String>
                }
                """);
            assertThat(TestUtils.getErrors(source)).isEmpty();
            assertThat(TestUtils.getWarnings(source)).isEmpty();
        }
    }

    @Nested
    class Callables {

        @Test
        void reportUndefinedFunctionCall() {
            assertErrors(check("""
                hello('world')
                """),
                tuple(1, 1, "`hello` is not defined"));
        }

        @Test
        void reportProcessOrWorkflowUsedAsVariableInTypedScript() {
            assertErrors(check("""
                nextflow.enable.types = true

                process hello {
                    exec:
                    true
                }

                workflow greet {
                }

                workflow {
                    channel.empty().mix(hello)
                    def x = greet
                }
                """),
                tuple(12, 25, "Process `hello` cannot be used as a variable"),
                tuple(13, 13, "Workflow `greet` cannot be used as a variable"));
        }

        @Test
        void allowProcessAsVariableInUntypedScript() {
            assertThat(check("""
                process hello {
                    exec:
                    true
                }

                workflow {
                    channel.empty().mix(hello)
                }
                """)).isEmpty();
        }

        @Test
        void reportProcessOutPropertyInTypedScript() {
            assertErrors(check("""
                nextflow.enable.types = true

                process hello {
                    output:
                    stdout()

                    script:
                    ''
                }

                workflow {
                    hello()
                    hello.out
                }
                """),
                tuple(13, 5, "Process `hello` cannot be used as a variable"));
        }

        @Test
        void reportProcessCalledFromFunction() {
            assertErrors(check("""
                process hello {
                    script:
                    \"""
                    echo hello
                    \"""
                }

                def greet() {
                    hello()
                }
                """),
                tuple(9, 5, "Processes can only be called from a workflow"));
        }

        @Test
        void reportProcessCalledFromClosure() {
            assertErrors(check("""
                process hello {
                    script:
                    \"""
                    echo hello
                    \"""
                }

                workflow {
                    [1, 2, 3].each { i ->
                        hello()
                    }
                }
                """),
                tuple(10, 9, "Processes cannot be called from within a closure"));
        }

        @Test
        void reportWorkflowCalledFromFunction() {
            assertErrors(check("""
                workflow hello {
                }

                def greet() {
                    hello()
                }
                """),
                tuple(5, 5, "Workflows can only be called from a workflow"));
        }

        @Test
        void reportWorkflowCalledFromClosure() {
            assertErrors(check("""
                workflow hello {
                }

                workflow {
                    [1, 2, 3].each { i ->
                        hello()
                    }
                }
                """),
                tuple(6, 9, "Workflows cannot be called from within a closure"));
        }

        @Test
        void resolvePipedOperatorsInEntryWorkflowWithParamsBlock() {
            assertThat(check("""
                params {
                    greeting: String = 'hello'
                }

                workflow {
                    channel.of(1, 2) | view
                }
                """)).isEmpty();
        }
    }

    @Nested
    class ReturnStatements {

        @Test
        void reportReturnInWorkflowBody() {
            assertErrors(check("""
                workflow {
                    main:
                    if( params.help )
                        return

                    publish:
                    nums = channel.of(1, 2, 3)
                }
                """),
                tuple(4, 9, "Return statement cannot be used in a workflow body -- use `exit()` instead"));
        }

        @Test
        void allowReturnInClosure() {
            assertThat(check("""
                workflow {
                    main:
                    nums = channel.of(1, 2, 3).map { v ->
                        if( v == 2 )
                            return 0
                        return v
                    }

                    publish:
                    out = nums
                }
                """)).isEmpty();
        }
    }

    @Nested
    class IncludedCallables {

        private static final String UTILS = "def basename(x) { 'hi' }\n";

        @TempDir
        Path root;

        @Test
        void doNotResolveIncludedFunctionAsVariable() {
            var main = TestUtils.tempFile(root, "main.nf", """
                include { basename } from './utils'

                process FOO {
                    shell:
                    "echo ${basename}"
                }

                workflow { FOO() }
                """);
            var utils = TestUtils.tempFile(root, "utils.nf", UTILS);
            assertThat(check(root, main, utils))
                .extracting(SyntaxException::getStartLine, SyntaxException::getOriginalMessage)
                .containsExactly(tuple(5, "`basename` is not defined"));
        }

        @Test
        void implicitDeclarationInProcessShadowsIncludedName() {
            var main = TestUtils.tempFile(root, "main.nf", """
                include { basename } from './utils'

                process FOO {
                    shell:
                    basename = basename(task)
                    "echo ${basename}"
                }

                workflow { FOO() }
                """);
            var utils = TestUtils.tempFile(root, "utils.nf", UTILS);
            assertThat(check(root, main, utils)).isEmpty();
        }

        @Test
        void explicitDeclarationShadowsIncludedFunction() {
            var main = TestUtils.tempFile(root, "main.nf", """
                include { basename } from './utils'

                process FOO {
                    shell:
                    def basename = basename(task)
                    "echo ${basename}"
                }

                workflow { FOO() }
                """);
            var utils = TestUtils.tempFile(root, "utils.nf", UTILS);
            assertThat(check(root, main, utils)).isEmpty();
        }

        @Test
        void resolveIncludedProcessAsVariableInWorkflow() {
            var main = TestUtils.tempFile(root, "main.nf", """
                include { BAR } from './modules/bar'

                workflow {
                    BAR()
                    BAR.out.view()
                }
                """);
            var module = TestUtils.tempFile(root, "modules/bar/main.nf", """
                process BAR {
                    output:
                    stdout

                    script:
                    "echo bar"
                }
                """);
            assertThat(check(root, main, module)).isEmpty();
        }
    }

}
