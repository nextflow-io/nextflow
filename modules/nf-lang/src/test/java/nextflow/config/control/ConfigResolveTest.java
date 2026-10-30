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
package nextflow.config.control;

import java.nio.file.Path;
import java.util.List;

import nextflow.test.TestUtils;
import org.codehaus.groovy.syntax.SyntaxException;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.tuple;

/**
 * Specifies the name resolution errors reported by {@link ConfigResolveVisitor}.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class ConfigResolveTest {

    private static List<SyntaxException> check(String contents) {
        return TestUtils.check(new ConfigParser(), contents);
    }

    private static List<SyntaxException> check(List<Path> files) {
        return TestUtils.check(new ConfigParser(), files);
    }

    @Test
    void reportUndefinedVariable() {
        var errors = check("""
            process.clusterOptions = "--cpus ${process.cpus}"
            """);
        assertThat(errors)
            .extracting(SyntaxException::getStartLine, SyntaxException::getStartColumn, SyntaxException::getOriginalMessage)
            .containsExactly(tuple(1, 36, "`process` is not defined"));
    }

    @Test
    void reportUndefinedVariableWithEnvHint() {
        var errors = check("""
            process.clusterOptions = "--cpus $PROCESS_CPUS"
            """);
        assertThat(errors)
            .extracting(SyntaxException::getStartLine, SyntaxException::getStartColumn, SyntaxException::getOriginalMessage)
            .containsExactly(tuple(1, 34, "`PROCESS_CPUS` is not defined (hint: use `env('...')` to access environment variable)"));
    }

    @Test
    void reportConfigIncludeInProfilesBlock(@TempDir Path root) {
        var config = TestUtils.tempFile(root, "nextflow.config", """
            profiles {
                includeConfig 'extra.config'
            }
            """);
        var extra = TestUtils.tempFile(root, "extra.config", "");
        var errors = check(List.of(config, extra));
        assertThat(errors)
            .extracting(SyntaxException::getStartLine, SyntaxException::getStartColumn, SyntaxException::getOriginalMessage)
            .containsExactly(tuple(2, 5, "Config includes are only allowed at the top-level or in a profile"));
    }

    @Test
    void allowConfigIncludeInProfile(@TempDir Path root) {
        var config = TestUtils.tempFile(root, "nextflow.config", """
            profiles {
                extra {
                    includeConfig 'extra.config'
                }
            }
            """);
        var extra = TestUtils.tempFile(root, "extra.config", "");
        assertThat(check(List.of(config, extra))).isEmpty();
    }

    @Test
    void reportUndeclaredVariableInClosure() {
        var errors = check("""
            process {
                clusterOptions = {
                    args_list = []
                    args_list.join(' ')
                }()
            }
            """);
        assertThat(errors)
            .extracting(SyntaxException::getStartLine, SyntaxException::getStartColumn, SyntaxException::getOriginalMessage)
            .containsExactly(
                tuple(3, 9, "`args_list` was assigned but not declared"),
                tuple(4, 9, "`args_list` is not defined"));
    }

    @Test
    void allowDeclaredVariableInClosure() {
        var errors = check("""
            process {
                clusterOptions = {
                    def args_list = []
                    args_list.join(' ')
                }()
            }
            """);
        assertThat(errors).isEmpty();
    }

    @Test
    void allowDynamicProcessDirectivesToReferenceProcessInputs() {
        var errors = check("""
            process {
                ext.prefix = { "${meta.id}.filter1" }

                publishDir = [
                    path: { "${params.outdir}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}" },
                    mode: params.publish_dir_mode,
                ]

                publishDir = [
                    path: { "${params.outdir}/imputation/${meta.tools}/samples/" },
                    mode: params.publish_dir_mode,
                ]
            }
            """);
        assertThat(errors).isEmpty();
    }

    @Test
    void reportProcessInputReferenceOutsideDynamicDirective() {
        var errors = check("""
            process {
                ext.prefix = { "${meta.id}.filter1" }()
            }
            """);
        assertThat(errors)
            .extracting(SyntaxException::getStartLine, SyntaxException::getStartColumn, SyntaxException::getOriginalMessage)
            .containsExactly(tuple(2, 23, "`meta` is not defined"));
    }

    @Test
    void allowTaskReferenceInExecutorJobName() {
        var errors = check("""
            executor.jobName = { "$task.name - $task.hash" }
            """);
        assertThat(errors).isEmpty();
    }

}
