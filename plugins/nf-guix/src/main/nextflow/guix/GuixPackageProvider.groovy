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

package nextflow.guix

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.packages.PackageProvider
import nextflow.packages.PackageSpec
import nextflow.util.Escape

/**
 * GNU Guix package provider implementation.
 *
 * Creates and activates GNU Guix-managed environments from a package list.
 */
@Slf4j
@CompileStatic
class GuixPackageProvider implements PackageProvider {

    private GuixCache cache
    private final GuixConfig config

    GuixPackageProvider(GuixConfig config) {
        this.config = config
        this.cache = new GuixCache(config)
    }

    @Override
    String getName() {
        return 'guix'
    }

    @Override
    boolean isAvailable() {
        try {
            def process = new ProcessBuilder('guix', '--version').start()
            process.waitFor()
            return process.exitValue() == 0
        } catch (Exception e) {
            log.debug "guix not available: ${e.message}"
            return false
        }
    }

    @Override
    Path createEnvironment(PackageSpec spec) {
        if (!supportsSpec(spec)) {
            throw new IllegalArgumentException("Unsupported package spec for guix: ${spec}")
        }

        String guixEnv
        if (spec.hasEnvironmentFile()) {
            guixEnv = spec.environment
        } else if (spec.hasEntries()) {
            // Space-separated package list e.g. 'bwa samtools'
            guixEnv = spec.entries.join(' ')
        } else {
            throw new IllegalArgumentException("Package spec must have either environment file or entries")
        }

        // per-process `options: [installOptions: '...']` overrides guix.installOptions
        final installOptionsOverride = spec.options?.get('installOptions') as String
        return cache.getCachePathFor(guixEnv, installOptionsOverride)
    }

    @Override
    String getActivationScript(Path envPath) {
        return """\
            export GUIX_PROFILE=${Escape.path(envPath)}
            . "\$GUIX_PROFILE/etc/profile"
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
        return ['manifest.scm', 'guix.scm']
    }
}
