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
import org.codehaus.groovy.syntax.SyntaxException;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.tuple;

/**
 * Specifies how {@link ResolveIncludeVisitor} resolves script includes.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class ResolveIncludeTest {

    @TempDir
    Path root;

    private List<SyntaxException> check(Path... files) {
        return TestUtils.check(new ScriptParser(root), List.of(files));
    }

    private Path file(String name, String contents) {
        return TestUtils.tempFile(root, name, contents);
    }

    @Test
    void reportInvalidIncludeSource() {
        var main = file("main.nf", """
            include { hello } from './hello.nf'
            """);
        var module = root.resolve("hello.nf");
        assertThat(check(main))
            .extracting(SyntaxException::getStartLine, SyntaxException::getStartColumn, SyntaxException::getOriginalMessage)
            .containsExactly(tuple(1, 1, "Invalid include source: '" + module + "'"));
    }

    @Test
    void reportIncludeSourceThatCouldNotBeParsed() {
        var main = file("main.nf", """
            include { hello } from './modules/hello'
            """);
        var module = file("modules/hello/main.nf", """
            println 1 2 3
            """);
        assertThat(check(main, module))
            .extracting(SyntaxException::getSourceLocator, SyntaxException::getStartLine, SyntaxException::getStartColumn, SyntaxException::getOriginalMessage)
            .containsExactly(
                tuple(main.toString(), 1, 1, "Module could not be parsed: '" + module + "'"),
                tuple(module.toString(), 1, 11, "Unexpected input: '2'"));
    }

    @Test
    void reportUndefinedInclude() {
        var main = file("main.nf", """
            include { hello } from './modules/hello'
            """);
        var module = file("modules/hello/main.nf", """
            workflow goodbye {
            }
            """);
        assertThat(check(main, module))
            .extracting(SyntaxException::getSourceLocator, SyntaxException::getStartLine, SyntaxException::getStartColumn, SyntaxException::getOriginalMessage)
            .containsExactly(tuple(main.toString(), 1, 1, "Included name 'hello' is not defined in module '" + module + "'"));
    }

    @Test
    void reportAliasedTypeInclude() {
        var main = file("main.nf", """
            include { Thing as Alias } from './types.nf'
            """);
        var module = file("types.nf", """
            record Thing {
                id: String
            }
            """);
        assertThat(check(main, module))
            .extracting(SyntaxException::getSourceLocator, SyntaxException::getStartLine, SyntaxException::getOriginalMessage)
            .containsExactly(tuple(main.toString(), 1, "Included types cannot be aliased"));
    }

    @ParameterizedTest
    @ValueSource(strings = { "workflow as GREET", "params as GreetParams" })
    void requireStaticTypingInScriptThatIncludesPipeline(String include) {
        var main = file("main.nf", """
            include { %s } from './greet.nf'
            """.formatted(include));
        var module = file("greet.nf", """
            nextflow.enable.types = true

            params {
                greeting: String = 'Hello'
            }

            workflow {
                println params.greeting
            }
            """);
        assertThat(check(main, module))
            .extracting(SyntaxException::getSourceLocator, SyntaxException::getOriginalMessage)
            .containsExactly(tuple(main.toString(), "Including a pipeline requires `nextflow.enable.types = true` in the including script"));
    }

    @Test
    void requireStaticTypingInIncludedPipeline() {
        var main = file("main.nf", """
            nextflow.enable.types = true

            include { workflow as GREET } from './greet.nf'
            """);
        var module = file("greet.nf", """
            params {
                greeting: String = 'Hello'
            }

            workflow {
                println params.greeting
            }
            """);
        assertThat(check(main, module))
            .extracting(SyntaxException::getSourceLocator, SyntaxException::getOriginalMessage)
            .containsExactly(tuple(main.toString(), "An included pipeline must enable static typing -- set `nextflow.enable.types = true` in '" + module + "'"));
    }

    @Test
    void resolveInclude() {
        var main = file("main.nf", """
            include { hello } from './modules/hello'
            """);
        var module = file("modules/hello/main.nf", """
            workflow hello {
            }
            """);
        assertThat(check(main, module)).isEmpty();
    }

    @Test
    void doNotCheckArgumentsOfPluginFunctionInTypedScript() {
        var main = file("main.nf", """
            nextflow.enable.types = true

            include { samplesheetToList } from 'plugin/nf-schema'

            workflow {
                def rows = samplesheetToList('s.csv', 'schema.json')
                println rows
                println rows.size()
            }
            """);
        assertThat(check(main)).isEmpty();
    }

    /**
     * An agent must bind through the ordinary include statement, in both the
     * plain-file and the directory ({@code main.nf}) form. The directory form
     * matters because getLocalIncludeUri is duplicated in ResolveIncludeVisitor
     * and ModuleResolver (lint/LSP vs run); only a test driving the include
     * catches a divergence.
     */
    @Nested
    class AgentModules {

        private static final String AGENT_MODULE = """
            nextflow.enable.types = true

            agent reporter {
                model 'openai/gpt-4o'
                instruction 'You write QA reports.'

                input:
                sample: String

                output:
                stdout()

                prompt:
                \"""
                Report on ${sample}.
                \"""
            }
            """;

        @Test
        void resolveAgentIncludeFromModuleFile() {
            var main = file("main.nf", """
                include { reporter } from './reporter.nf'

                workflow {
                    reporter('s1')
                }
                """);
            var module = file("reporter.nf", AGENT_MODULE);
            assertThat(check(main, module)).isEmpty();
        }

        @Test
        void resolveAgentIncludeFromModuleDirectory() {
            var main = file("main.nf", """
                include { reporter } from './mods/reporter'

                workflow {
                    reporter('s1')
                }
                """);
            var module = file("mods/reporter/main.nf", AGENT_MODULE);
            assertThat(check(main, module)).isEmpty();
        }

        @Test
        void resolveAliasedAgentIncludeAtCallSite() {
            var main = file("main.nf", """
                include { reporter as qc } from './mods/reporter'

                workflow {
                    qc('s1')
                }
                """);
            var module = file("mods/reporter/main.nf", AGENT_MODULE);
            assertThat(check(main, module)).isEmpty();
        }

        @Test
        void resolveAgentIncludeThatIsNeverCalled() {
            var main = file("main.nf", """
                include { reporter } from './mods/reporter'
                """);
            var module = file("mods/reporter/main.nf", AGENT_MODULE);
            assertThat(check(main, module)).isEmpty();
        }

        @Test
        void reportUndefinedAgentInclude() {
            var main = file("main.nf", """
                include { nope } from './mods/reporter'
                """);
            var module = file("mods/reporter/main.nf", AGENT_MODULE);
            assertThat(check(main, module))
                .extracting(SyntaxException::getSourceLocator, SyntaxException::getOriginalMessage)
                .containsExactly(tuple(main.toString(), "Included name 'nope' is not defined in module '" + module + "'"));
        }
    }

}
