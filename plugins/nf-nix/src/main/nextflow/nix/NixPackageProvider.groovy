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

package nextflow.nix

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.packages.PackageProvider
import nextflow.packages.PackageSpec
import nextflow.util.Escape

/**
 * Nix package provider implementation.
 *
 * Creates and activates Nix-managed environments from a package list
 * or an existing Nix profile directory.
 */
@Slf4j
@CompileStatic
class NixPackageProvider implements PackageProvider {

    private NixCache cache
    private final NixConfig config

    NixPackageProvider(NixConfig config) {
        this.config = config
        this.cache = new NixCache(config)
    }

    @Override
    String getName() {
        return 'nix'
    }

    @Override
    boolean isAvailable() {
        try {
            def process = new ProcessBuilder('nix', '--version').start()
            process.waitFor()
            return process.exitValue() == 0
        } catch (Exception e) {
            log.debug "nix not available: ${e.message}"
            return false
        }
    }

    @Override
    Path createEnvironment(PackageSpec spec) {
        if (!supportsSpec(spec)) {
            throw new IllegalArgumentException("Unsupported package spec for nix: ${spec}")
        }

        String env
        if (spec.hasEnvironmentFile()) {
            // Path to an existing Nix profile directory
            env = spec.environment
        } else if (spec.hasEntries()) {
            // Space-separated package list e.g. 'bwa samtools'
            env = spec.entries.join(' ')
        } else {
            throw new IllegalArgumentException("Package spec must have either environment file or entries")
        }

        // per-process `options: [installOptions: '...']` overrides nix.installOptions
        final installOptionsOverride = spec.options?.get('installOptions') as String
        return cache.getCachePathFor(env, installOptionsOverride)
    }

    @Override
    String getActivationScript(Path envPath) {
        return """\
            export PATH=${Escape.path(envPath)}/bin:\$PATH
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
}
