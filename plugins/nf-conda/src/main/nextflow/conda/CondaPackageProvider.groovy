/*
 * Copyright 2013-2024, Seqera Labs
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

package nextflow.conda

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.packages.PackageProvider
import nextflow.packages.PackageSpec
import nextflow.util.Escape

/**
 * Conda package provider implementation
 * 
 * @author Edmund Miller <edmund.miller@seqera.io>
 */
@Slf4j
@CompileStatic
class CondaPackageProvider implements PackageProvider {

    private final CondaCache cache
    private final CondaConfig config

    CondaPackageProvider(CondaConfig config) {
        this.config = config
        this.cache = new CondaCache(config)
    }

    @Override
    String getName() {
        return 'conda'
    }

    @Override
    boolean isAvailable() {
        try {
            def process = new ProcessBuilder('conda', '--version').start()
            process.waitFor()
            return process.exitValue() == 0
        } catch (Exception e) {
            log.debug "Conda not available: ${e.message}"
            return false
        }
    }

    @Override
    Path createEnvironment(PackageSpec spec) {
        if (!supportsSpec(spec)) {
            throw new IllegalArgumentException("Unsupported package spec for conda: ${spec}")
        }

        String condaEnv
        if (spec.hasEnvironmentFile()) {
            // Handle environment file
            condaEnv = spec.environment
        } else if (spec.hasEntries()) {
            // Handle package list
            condaEnv = spec.entries.join(' ')
        } else {
            throw new IllegalArgumentException("Package spec must have either environment file or entries")
        }

        return cache.getCachePathFor(condaEnv)
    }

    @Override
    String getActivationScript(Path envPath) {
        def binaryName = cache.getBinaryName()
        final command = config.useMicromamba()
            ? 'eval "$(micromamba shell hook --shell bash)" && micromamba activate'
            : 'source $(conda info --json | awk \'/conda_prefix/ { gsub(/"|,/, "", $2); print $2 }\')/bin/activate'
        
        return """\
            ${command} ${Escape.path(envPath)}
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