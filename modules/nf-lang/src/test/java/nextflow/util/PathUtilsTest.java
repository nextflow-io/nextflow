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
package nextflow.util;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.HashSet;
import java.util.List;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;

import static nextflow.test.TestUtils.tempFile;
import static org.assertj.core.api.Assertions.assertThat;

/**
 * Specifies the file visiting and exclusion helpers.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class PathUtilsTest {

    @Test
    void visitAllFilesThatSatisfyCondition(@TempDir Path root) throws IOException {
        var files = List.of(
            "main.nf",
            "nextflow.config",
            "modules/index/main.nf",
            "modules/fastqc/main.nf",
            "modules/quant/main.nf",
            "modules/multiqc/main.nf",
            "modules/rnaseq.nf"
        );
        for( var name : files )
            tempFile(root, name);

        var result = new HashSet<String>();
        PathUtils.visitFiles(
            root,
            path -> Files.isDirectory(path) || path.toString().endsWith(".nf"),
            path -> result.add(root.relativize(path).toString()));

        assertThat(result).containsExactlyInAnyOrder(
            "main.nf",
            "modules/index/main.nf",
            "modules/fastqc/main.nf",
            "modules/quant/main.nf",
            "modules/multiqc/main.nf",
            "modules/rnaseq.nf"
        );
    }

    @ParameterizedTest
    @CsvSource(delimiter = '|', textBlock = """
        .git                        | true
        .nf-test                    | true
        modules/foo/.nf-test        | true
        modules/foo/bar             | true
        work                        | true
        nf-test.config              | true
        modules/foo/nf-test.config  | true
        main.nf                     | false
        modules/foo/main.nf         | false
        subworkflows/foobar/main.nf | false
        nextflow.config             | false
        conf/base.config            | false
        """)
    void determineWhetherPathIsExcluded(String path, boolean excluded) {
        var patterns = List.of(".git", ".nf-test", "work", "foo/bar", "nf-test.config");
        assertThat(PathUtils.isExcluded(Path.of(path), patterns)).isEqualTo(excluded);
    }

}
