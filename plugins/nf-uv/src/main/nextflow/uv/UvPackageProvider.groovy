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

package nextflow.uv

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.packages.PackageProvider
import nextflow.packages.PackageSpec
import nextflow.util.Escape

/**
 * uv package provider implementation.
 *
 * Creates and activates uv-managed Python virtual environments from a
 * package list, a requirements file, a pyproject.toml file, or an existing
 * virtual environment directory.
 */
@Slf4j
@CompileStatic
class UvPackageProvider implements PackageProvider {

    private UvCache cache
    private final UvConfig config

    UvPackageProvider(UvConfig config) {
        this.config = config
        this.cache = new UvCache(config)
    }

    @Override
    String getName() {
        return 'uv'
    }

    @Override
    boolean isAvailable() {
        try {
            def process = new ProcessBuilder('uv', '--version').start()
            process.waitFor()
            return process.exitValue() == 0
        } catch (Exception e) {
            log.debug "uv not available: ${e.message}"
            return false
        }
    }

    @Override
    Path createEnvironment(PackageSpec spec) {
        if (!supportsSpec(spec)) {
            throw new IllegalArgumentException("Unsupported package spec for uv: ${spec}")
        }

        String uvEnv
        if (spec.hasEnvironmentFile()) {
            // Path to a requirements file, pyproject.toml, or existing venv directory
            uvEnv = spec.environment
        } else if (spec.hasEntries()) {
            // Space-separated package list e.g. 'numpy pandas matplotlib',
            // or a single path to a requirements/pyproject file or venv directory
            uvEnv = spec.entries.join(' ')
        } else {
            throw new IllegalArgumentException("Package spec must have either environment file or entries")
        }

        // per-process `options: [installOptions: '...']` overrides uv.installOptions
        final installOptionsOverride = spec.options?.get('installOptions') as String
        return cache.getCachePathFor(uvEnv, installOptionsOverride)
    }

    @Override
    String getActivationScript(Path envPath) {
        return """\
            source ${Escape.path(envPath)}/bin/activate
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
        return ['requirements.txt', 'pyproject.toml']
    }
}
