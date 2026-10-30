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
 * Specifies how {@link ResolveIncludeVisitor} resolves config includes.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class ResolveIncludeTest {

    @TempDir
    Path root;

    private static List<SyntaxException> check(List<Path> files) {
        return TestUtils.check(new ConfigParser(), files);
    }

    @Test
    void reportInvalidIncludeSource() {
        var config = TestUtils.tempFile(root, "nextflow.config", """
            includeConfig 'extra.config'
            """);
        var extra = root.resolve("extra.config");
        assertThat(check(List.of(config)))
            .extracting(SyntaxException::getStartLine, SyntaxException::getStartColumn, SyntaxException::getOriginalMessage)
            .containsExactly(tuple(1, 1, "Invalid include source: '" + extra + "'"));
    }

    @Test
    void resolveInclude() {
        var config = TestUtils.tempFile(root, "nextflow.config", """
            includeConfig 'extra.config'
            """);
        var extra = TestUtils.tempFile(root, "extra.config", """
            process.cpus = 8
            """);
        assertThat(check(List.of(config, extra))).isEmpty();
    }

    @Test
    void doNotResolveRemoteInclude() {
        var config = TestUtils.tempFile(root, "nextflow.config", """
            includeConfig 'http://example.com/nextflow.config'
            """);
        assertThat(check(List.of(config))).isEmpty();
    }

    @Test
    void doNotResolveDynamicInclude() {
        var config = TestUtils.tempFile(root, "nextflow.config", """
            includeConfig ({ 'extra.config' }())
            """);
        assertThat(check(List.of(config))).isEmpty();
    }

}
