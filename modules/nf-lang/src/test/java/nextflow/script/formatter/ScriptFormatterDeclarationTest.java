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
package nextflow.script.formatter;

import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;

import static nextflow.script.formatter.FormatterTestUtils.DEFAULT;
import static nextflow.script.formatter.FormatterTestUtils.assertFormat;
import static nextflow.script.formatter.FormatterTestUtils.assertFormatted;
import static nextflow.script.formatter.FormatterTestUtils.format;
import static org.assertj.core.api.Assertions.assertThat;

/**
 * Specifies the output of {@link ScriptFormattingVisitor} for script
 * declarations.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class ScriptFormatterDeclarationTest {

    @Test
    void formatCodeSnippet() {
        assertFormat(
            """
            println 'Hello!'
            """,
            """
            println('Hello!')
            """
        );
    }

    @Test
    void doNotSortScriptDeclarationsByDefault() {
        assertFormatted("""
            params.foo = 'bar'

            include { foo ; bar } from './foobar.nf'

            process hello {
                script:
                'echo true'
            }

            workflow hello {
            }

            workflow {
            }

            output {
            }
            """);
    }

    @Nested
    class Includes {

        @Test
        void formatSingleLineInclude() {
            assertFormat(
                """
                include{foo;bar}from'./foobar.nf'
                """,
                """
                include { foo ; bar } from './foobar.nf'
                """
            );
        }

        @Test
        void formatMultiLineInclude() {
            assertFormat(
                """
                include{
                foo;bar
                }from'./foobar.nf'
                """,
                """
                include {
                    foo ;
                    bar
                } from './foobar.nf'
                """
            );
        }
    }

    @Nested
    class Params {

        @Test
        void formatParameterDeclaration() {
            assertFormat(
                """
                params.foo='bar'
                """,
                """
                params.foo = 'bar'
                """
            );
        }
    }

    @Nested
    class Workflows {

        @Test
        void formatLegacyWorkflow() {
            assertFormat(
                """
                workflow hello{
                take: x ; y ; main: xy=x*y ; emit: result = xy
                }
                """,
                """
                workflow hello {
                    take:
                    x
                    y

                    main:
                    xy = x * y

                    emit:
                    result = xy
                }
                """
            );
        }

        @Test
        void formatLegacyWorkflowWithoutMain() {
            assertFormat(
                """
                workflow hello{
                take: x ; y ; emit: result = x * y
                }
                """,
                """
                workflow hello {
                    take:
                    x
                    y

                    emit:
                    result = x * y
                }
                """
            );
        }

        @Test
        void formatTypedWorkflow() {
            assertFormat(
                """
                nextflow.enable.types = true

                workflow hello{
                take: x:Integer ; y:Integer ; main: xy=x*y ; emit: result:Integer = xy
                }
                """,
                """
                nextflow.enable.types = true

                workflow hello {
                    take:
                    x: Integer
                    y: Integer

                    main:
                    xy = x * y

                    emit:
                    result: Integer = xy
                }
                """
            );
        }

        @Test
        void doNotCopyOutputTypeIntoPublishTarget() {
            var result = format("""
                nextflow.enable.types = true

                workflow {
                    main:
                    ch = channel.of(1)

                    publish:
                    a = ch
                    b: Channel<String> = ch
                }

                output {
                    a: Channel<Integer> {
                        path 'a'
                    }
                    b: Channel<Integer> {
                        path 'b'
                    }
                }
                """, DEFAULT, true);
            assertThat(result).contains("    a = ch\n", "    b = ch\n");
        }
    }

    @Nested
    class Processes {

        @Test
        void formatLegacyProcess() {
            assertFormat(
                """
                process hello{
                debug(true) ; input: val x ; val y ; output: tuple val(x), val(y) ; script: 'echo true'
                }
                """,
                """
                process hello {
                    debug true

                    input:
                    val x
                    val y

                    output:
                    tuple val(x), val(y)

                    script:
                    'echo true'
                }
                """
            );
        }

        @Test
        void formatTypedProcess() {
            assertFormat(
                """
                nextflow.enable.types=true

                process hello{
                debug(true) ; input: tuple(id:String,infile:Path) ; index:Path ; stage: stageAs(infile,'input.txt') ; output: result=tuple(id,file('output.txt')) ; script: 'cat input.txt > output.txt'
                }
                """,
                """
                nextflow.enable.types = true

                process hello {
                    debug true

                    input:
                    tuple(id: String, infile: Path)
                    index: Path

                    stage:
                    stageAs infile, 'input.txt'

                    output:
                    result = tuple(id, file('output.txt'))

                    script:
                    'cat input.txt > output.txt'
                }
                """
            );
        }

        @Test
        void formatRecordInputOnSeparateLines() {
            assertFormat(
                """
                nextflow.enable.types=true

                process hello{
                input: record(id:String,infile:Path) ; script: 'cat input.txt > output.txt'
                }
                """,
                """
                nextflow.enable.types = true

                process hello {
                    input:
                    record(
                        id: String,
                        infile: Path
                    )

                    script:
                    'cat input.txt > output.txt'
                }
                """
            );
        }
    }

    @Nested
    class Agents {

        // an agent declaration used to be walked without emitting anything, so
        // formatting an agent module deleted its entire content
        @Test
        void formatAgent() {
            assertFormatted("""
                nextflow.enable.types = true

                agent reporter {
                    model 'openai/gpt-4o'
                    instruction 'You write QA reports.'
                    skills 'qa-report', 'style'

                    input:
                    sample: String
                    depth: Integer

                    output:
                    stdout()

                    prompt:
                    \"""
                    Write a report for ${sample} at depth ${depth}.
                    \"""
                }
                """);
        }
    }

    @Nested
    class Functions {

        @Test
        void formatUntypedFunction() {
            assertFormat(
                """
                def hello(x,y){
                def xy=x*y ; return xy
                }
                """,
                """
                def hello(x, y) {
                    def xy = x * y
                    return xy
                }
                """
            );
        }

        @Test
        void convertLegacyTypeAnnotations() {
            assertFormat(
                """
                Integer hello(Integer x,Integer y){
                Integer xy=x*y ; return xy
                }
                """,
                """
                def hello(x: Integer, y: Integer) -> Integer {
                    def xy: Integer = x * y
                    return xy
                }
                """
            );
        }

        @Test
        void formatTypedFunction() {
            assertFormat(
                """
                def hello(x:Integer,y:Integer)->Integer{
                def xy=x*y ; return xy
                }
                """,
                """
                def hello(x: Integer, y: Integer) -> Integer {
                    def xy = x * y
                    return xy
                }
                """
            );
        }

        @Test
        void dropObjectAndArrayTypes() {
            assertFormat(
                """
                Object hello(Object x,String[] y){
                Object xy=x ; return xy
                }
                """,
                """
                def hello(x, y) {
                    def xy = x
                    return xy
                }
                """
            );
        }
    }

    @Nested
    class EnumsAndRecords {

        @Test
        void formatEnum() {
            assertFormat(
                """
                enum Colors{RED,GREEN,BLUE}
                """,
                """
                enum Colors {
                    RED,
                    GREEN,
                    BLUE,
                }
                """
            );
        }

        @Test
        void formatRecord() {
            assertFormat(
                """
                record FastqPair{id:String;fastq_1: Path;fastq_2: Path?}
                """,
                """
                record FastqPair {
                    id: String
                    fastq_1: Path
                    fastq_2: Path?
                }
                """
            );
        }
    }

    @Nested
    class Outputs {

        @Test
        void formatOutputBlock() {
            assertFormat(
                """
                workflow{}

                output{
                foo{path'foo'}
                bar{path'bar';index{path'index.csv'}}
                }
                """,
                """
                workflow {
                }

                output {
                    foo {
                        path 'foo'
                    }
                    bar {
                        path 'bar'
                        index {
                            path 'index.csv'
                        }
                    }
                }
                """
            );
        }

        @Test
        void formatTypedOutputs() {
            assertFormat(
                """
                workflow{}

                output{
                foo:Path{path'foo'}
                bar:Channel<Path>{path'bar';index{path'index.csv'}}
                }
                """,
                """
                workflow {
                }

                output {
                    foo: Path {
                        path 'foo'
                    }
                    bar: Channel<Path> {
                        path 'bar'
                        index {
                            path 'index.csv'
                        }
                    }
                }
                """
            );
        }

        @Test
        void formatIndexDirectives() {
            assertFormat(
                """
                workflow{}

                output{
                foo{index{header ( ['a','b'] );sep('\\t')}}
                }
                """,
                """
                workflow {
                }

                output {
                    foo {
                        index {
                            header(['a', 'b'])
                            sep '\\t'
                        }
                    }
                }
                """
            );
        }
    }

}
