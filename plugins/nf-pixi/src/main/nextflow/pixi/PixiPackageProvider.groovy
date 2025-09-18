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

package nextflow.pixi

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.packages.PackageProvider
import nextflow.packages.PackageSpec
import nextflow.util.Escape

/**
 * Pixi package provider implementation
 * 
 * @author Edmund Miller <edmund.miller@seqera.io>
 */
@Slf4j
@CompileStatic
class PixiPackageProvider implements PackageProvider {

    private final PixiCache cache
    private final PixiConfig config

    PixiPackageProvider(PixiConfig config) {
        this.config = config
        this.cache = new PixiCache(config)
    }

    @Override
    String getName() {
        return 'pixi'
    }

    @Override
    boolean isAvailable() {
        try {
            def process = new ProcessBuilder('pixi', '--version').start()
            process.waitFor()
            return process.exitValue() == 0
        } catch (Exception e) {
            log.debug "Pixi not available: ${e.message}"
            return false
        }
    }

    @Override
    Path createEnvironment(PackageSpec spec) {
        if (!supportsSpec(spec)) {
            throw new IllegalArgumentException("Unsupported package spec for pixi: ${spec}")
        }

        String pixiEnv
        if (spec.hasEnvironmentFile()) {
            // Handle environment file
            pixiEnv = spec.environment
        } else if (spec.hasEntries()) {
            // Handle package list
            pixiEnv = spec.entries.join(' ')
        } else {
            throw new IllegalArgumentException("Package spec must have either environment file or entries")
        }

        return cache.getCachePathFor(pixiEnv)
    }

    @Override
    String getActivationScript(Path envPath) {
        def result = ""

        // Check if there's a .pixi file that points to the project directory
        final pixiFile = envPath.resolve('.pixi')
        if (pixiFile.exists()) {
            // Read the project directory path
            final projectDir = pixiFile.text.trim()
            result += "cd ${Escape.path(projectDir as String)} && "
            result += "eval \"\$(pixi shell-hook --shell bash)\" && "
            result += "cd \"\$OLDPWD\"\n"
        } else {
            // Direct activation from environment directory
            result += "cd ${Escape.path(envPath)} && "
            result += "eval \"\$(pixi shell-hook --shell bash)\" && "
            result += "cd \"\$OLDPWD\"\n"
        }

        return result
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