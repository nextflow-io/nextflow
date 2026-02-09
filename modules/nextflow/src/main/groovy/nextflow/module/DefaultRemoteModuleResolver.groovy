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

package nextflow.module

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Global
import nextflow.NF
import nextflow.Session
import nextflow.config.ModulesConfig
import nextflow.config.RegistryConfig
import nextflow.exception.IllegalModulePath
import nextflow.module.spi.RemoteModuleResolver
import nextflow.pipeline.PipelineSpec

import java.nio.file.Path

/**
 * Default implementation of RemoteModuleResolver using the Nextflow module registry.
 *
 * <p>This implementation:
 * <ul>
 *   <li>Checks for locally installed modules in the project's modules directory</li>
 *   <li>Downloads modules from the configured registry if not present</li>
 *   <li>Reads version constraints from nextflow_spec.json</li>
 *   <li>Uses the Session's registry configuration</li>
 * </ul>
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Slf4j
@CompileStatic
class DefaultRemoteModuleResolver implements RemoteModuleResolver {

    @Override
    Path resolve(String moduleName, Path baseDir) {

        final modulesConfig = getModuleConfig(baseDir)
        final registryConfig = Global.config.navigate('registry') as RegistryConfig

        // Create module resolver
        def resolver = new ModuleResolver(baseDir, modulesConfig, registryConfig)

        try {
            log.debug "Resolving remote module: ${moduleName}"

            // Parse module reference
            def reference = ModuleReference.parse(moduleName)

            // Resolve module (will auto-install if missing or version mismatch)
            def mainFile = resolver.resolve(reference, null, true)

            log.info "Module ${reference.nameWithoutPrefix} resolved to ${mainFile}"
            return mainFile
        } catch (Exception e) {
            throw new IllegalModulePath(
                "Failed to resolve remote module ${moduleName}: ${e.message}",
                e
            )
        }
    }

    @Override
    int getPriority() {
        return 0  // Default implementation has lowest priority
    }

    private ModulesConfig getModuleConfig(Path baseDir) {
        def specFile = new PipelineSpec(baseDir)

        if (!specFile.exists()) {
            log.warn1("Remote module specified and 'nextflow_spec.json' not found.")
            return new ModulesConfig()
        }

        def modules = specFile.getModules()
        if (!modules || modules.isEmpty()) {
            log.warn1("Remote module specified and no modules configured in 'nextflow_spec.json'")
            return new ModulesConfig()
        }
        return new ModulesConfig(modules)
    }
}