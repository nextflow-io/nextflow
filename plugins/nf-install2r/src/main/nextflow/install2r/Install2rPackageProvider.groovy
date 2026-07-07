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

package nextflow.install2r

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.packages.PackageProvider
import nextflow.packages.PackageSpec
import nextflow.util.Escape

/**
 * install2.r package provider implementation.
 *
 * Installs R/CRAN packages into a library directory using install2.r (from the
 * littler package) and activates it by pointing R_LIBS_USER at that directory.
 */
@Slf4j
@CompileStatic
class Install2rPackageProvider implements PackageProvider {

    private Install2rCache cache
    private final Install2rConfig config

    Install2rPackageProvider(Install2rConfig config) {
        this.config = config
        this.cache = new Install2rCache(config)
    }

    @Override
    String getName() {
        return 'install2r'
    }

    @Override
    boolean isAvailable() {
        try {
            def process = new ProcessBuilder('bash', '-c', 'command -v install2.r >/dev/null 2>&1').start()
            process.waitFor()
            return process.exitValue() == 0
        } catch (Exception e) {
            log.debug "install2.r not available: ${e.message}"
            return false
        }
    }

    @Override
    Path createEnvironment(PackageSpec spec) {
        if (!supportsSpec(spec)) {
            throw new IllegalArgumentException("Unsupported package spec for install2.r: ${spec}")
        }

        String env
        if (spec.hasEnvironmentFile()) {
            env = spec.environment
        } else if (spec.hasEntries()) {
            env = spec.entries.join(' ')
        } else {
            throw new IllegalArgumentException("Package spec must have either environment file or entries")
        }

        // per-process `options: [installOptions: '...']` overrides install2r.installOptions
        final installOptionsOverride = spec.options?.get('installOptions') as String
        return cache.getCachePathFor(env, installOptionsOverride)
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
}
