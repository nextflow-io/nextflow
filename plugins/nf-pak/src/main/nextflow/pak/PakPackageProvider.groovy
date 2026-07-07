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

package nextflow.pak

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.packages.PackageProvider
import nextflow.packages.PackageSpec
import nextflow.util.Escape

/**
 * R pak package provider implementation.
 *
 * Creates and activates pak-managed R library directories from a package
 * list or an existing R library directory. Packages are installed with
 * {@code pak::pkg_install(..., lib=<dir>)} and activated by pointing
 * {@code R_LIBS_USER} at the library directory.
 */
@Slf4j
@CompileStatic
class PakPackageProvider implements PackageProvider {

    private PakCache cache
    private final PakConfig config

    PakPackageProvider(PakConfig config) {
        this.config = config
        this.cache = new PakCache(config)
    }

    @Override
    String getName() {
        return 'pak'
    }

    @Override
    boolean isAvailable() {
        try {
            def process = new ProcessBuilder('Rscript', '-e', 'if(!requireNamespace("pak",quietly=TRUE)) quit(status=1)').start()
            process.waitFor()
            return process.exitValue() == 0
        } catch (Exception e) {
            log.debug "R pak not available: ${e.message}"
            return false
        }
    }

    @Override
    Path createEnvironment(PackageSpec spec) {
        if (!supportsSpec(spec)) {
            throw new IllegalArgumentException("Unsupported package spec for pak: ${spec}")
        }

        String pakEnv
        if (spec.hasEnvironmentFile()) {
            // Path to an existing R library directory
            pakEnv = spec.environment
        } else if (spec.hasEntries()) {
            // Space-separated package list e.g. 'dplyr ggplot2',
            // or a single path to an existing R library directory
            pakEnv = spec.entries.join(' ')
        } else {
            throw new IllegalArgumentException("Package spec must have either environment file or entries")
        }

        // per-process `options: [installOptions: '...']` overrides pak.installOptions
        final installOptionsOverride = spec.options?.get('installOptions') as String
        return cache.getCachePathFor(pakEnv, installOptionsOverride)
    }

    @Override
    String getActivationScript(Path envPath) {
        return """\
            export R_LIBS_USER=${Escape.path(envPath)}
            """.stripIndent()
    }

    @Override
    boolean supportsSpec(PackageSpec spec) {
        return spec.provider == getName()
    }

    @Override
    Object getConfig() {
        return config
    }

    @Override
    List<String> getManifestFileNames() {
        // Only DESCRIPTION is supported via pak::local_install_deps(); renv.lock
        // is intentionally excluded because pak::lockfile_install() reads pak's
        // own native lock format, not the renv.lock JSON schema.
        return ['DESCRIPTION']
    }
}
