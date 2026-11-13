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
package nextflow.config.formatter;

import nextflow.config.control.ConfigParser;
import nextflow.script.formatter.FormattingOptions;
import org.junit.jupiter.api.Test;

import static org.assertj.core.api.Assertions.assertThat;

/**
 * Specifies the output of {@link ConfigFormattingVisitor}.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class ConfigFormatterTest {

    private static String format(String contents) {
        var source = new ConfigParser().parse("nextflow.config", contents);
        assertThat(source.getErrorCollector().getErrors()).isNullOrEmpty();
        var formatter = new ConfigFormattingVisitor(source, new FormattingOptions(4, true));
        formatter.visit();
        assertThat(formatter.getMissingComments()).isEmpty();
        return formatter.toString();
    }

    /**
     * Assert that {@code input} is formatted as {@code expected}, and that
     * {@code expected} is already formatted.
     */
    private static void assertFormat(String input, String expected) {
        assertThat(format(input)).isEqualTo(expected);
        assertThat(format(expected)).as("formatting is idempotent").isEqualTo(expected);
    }

    private static void assertFormatted(String source) {
        assertThat(format(source)).isEqualTo(source);
    }

    @Test
    void formatAssignment() {
        assertFormat(
            """
            process . clusterOptions={"--cpus ${task.cpus}"}
            """,
            """
            process.clusterOptions = { "--cpus ${task.cpus}" }
            """
        );
    }

    @Test
    void formatBlock() {
        assertFormat(
            """
            process{clusterOptions={"--cpus ${task.cpus}"}}
            """,
            """
            process {
                clusterOptions = { "--cpus ${task.cpus}" }
            }
            """
        );
    }

    @Test
    void formatApplyBlock() {
        assertFormat(
            """
            plugins{id'nf-hello'}
            """,
            """
            plugins {
                id 'nf-hello'
            }
            """
        );
    }

    @Test
    void preserveComments() {
        assertFormatted("""
            // about the process scope
            process {
                cpus = 2 // how many
                // dangling
            }

            // about the plugins
            plugins {
                id 'nf-hello' // a plugin
                // dangling
            }

            // trailing comment at the end of the file
            """);
    }

    @Test
    void preserveCommentOnInclude() {
        assertFormatted("""
            // about the include
            includeConfig 'other.config' // over here
            """);
    }

    @Test
    void doNotTreatCommentLikeStringTextAsComment() {
        assertFormatted("""
            params.pattern = "${params.dir}/*.vcf.gz"
            params.url = "${params.host}//path"
            """);
    }

    @Test
    void reindentMultiLineStringWithItsStatement() {
        assertFormat(
            """
            process {
              beforeScript = \"""
              echo one
              \"""
            }
            """,
            """
            process {
                beforeScript = \"""
                echo one
                \"""
            }
            """
        );
    }

}
