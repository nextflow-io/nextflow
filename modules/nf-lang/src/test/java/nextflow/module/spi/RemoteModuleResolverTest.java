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
package nextflow.module.spi;

import java.io.IOException;
import java.net.URI;
import java.nio.file.Files;
import java.nio.file.Path;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import static nextflow.module.spi.RemoteModuleResolver.resolveBaseDir;
import static org.assertj.core.api.Assertions.assertThat;

/**
 * Specifies the base directory rules of remote module include resolution.
 *
 * @author Ben Sherman <ben.sherman@seqera.io>
 */
class RemoteModuleResolverTest {

    @TempDir
    Path projectDir;

    @Test
    void resolveIncludeInModuleScriptAgainstModuleDirectory() {
        var moduleDir = projectDir.resolve("modules/nf-core/my-workflow");
        assertThat(resolveBaseDir(moduleDir.resolve("main.nf"), projectDir)).isEqualTo(moduleDir);

        // at any depth
        var nestedDir = moduleDir.resolve("modules/nf-core/star/align");
        assertThat(resolveBaseDir(nestedDir.resolve("main.nf"), projectDir)).isEqualTo(nestedDir);
    }

    @Test
    void resolveIncludeInNonModuleScriptAgainstProjectDirectory() {
        // an entry script in a subdirectory, as with a nested manifest.mainScript
        var scriptDir = projectDir.resolve("workflows");
        assertThat(resolveBaseDir(scriptDir.resolve("main.nf"), projectDir)).isEqualTo(projectDir);

        assertThat(resolveBaseDir(projectDir.resolve("main.nf"), projectDir)).isEqualTo(projectDir);
    }

    @Test
    void fallBackToProjectDirectoryWhenIncludingScriptIsUnknown() {
        assertThat(resolveBaseDir((Path) null, projectDir)).isEqualTo(projectDir);
        assertThat(resolveBaseDir((URI) null, projectDir)).isEqualTo(projectDir);
    }

    @Test
    void doNotTreatScriptOutsideProjectAsVendoredModule() {
        assertThat(resolveBaseDir(Path.of("/elsewhere/modules/a/b/main.nf"), projectDir)).isEqualTo(projectDir);
    }

    @Test
    void useIncludingScriptDirectoryWhenProjectDirectoryIsUnknown() {
        // as when generating a spec for a module that has no spec yet
        var moduleDir = projectDir.resolve("my-workflow");
        assertThat(resolveBaseDir(moduleDir.resolve("main.nf"), null)).isEqualTo(moduleDir);
    }

    @Test
    void resolveScriptIdentifiedByUri() throws IOException {
        var moduleDir = Files.createDirectories(projectDir.resolve("modules/nf-core/my-workflow"));
        assertThat(resolveBaseDir(moduleDir.resolve("main.nf").toUri(), projectDir)).isEqualTo(moduleDir);
    }

    @Test
    void resolveNonFileUriAgainstProjectDirectory() {
        assertThat(resolveBaseDir(URI.create("untitled:Untitled-1"), projectDir)).isEqualTo(projectDir);
    }

}
