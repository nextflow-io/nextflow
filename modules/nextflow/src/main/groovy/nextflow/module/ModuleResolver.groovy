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
import nextflow.config.ModulesConfig
import nextflow.config.RegistryConfig
import nextflow.exception.AbortOperationException

import java.nio.file.Files
import java.nio.file.Path

/**
 * Core module resolution logic that coordinates registry, storage, and version management
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Slf4j
@CompileStatic
class ModuleResolver {

    private final ModuleRegistryClient registryClient
    private final ModuleStorage storage
    private final ModulesConfig modulesConfig
    private final RegistryConfig registryConfig

    ModuleResolver(Path baseDir, ModulesConfig modulesConfig = null, RegistryConfig registryConfig = null) {
        this.registryConfig = registryConfig ?: new RegistryConfig()
        this.registryClient = new ModuleRegistryClient(this.registryConfig)
        this.storage = new ModuleStorage(baseDir)
        this.modulesConfig = modulesConfig ?: new ModulesConfig()
    }

    /**
     * Resolve a module reference to an installed module path
     *
     * @param reference The module reference
     * @param version Optional specific version (null = use config or latest)
     * @param autoInstall Whether to auto-install if not present (default: false)
     * @return Path to the module's main.nf file
     */
    Path resolve(ModuleReference reference, String version = null, boolean autoInstall = false) {
        // Determine version: explicit > config > latest
        def targetVersion = version ?: modulesConfig.getVersion(reference.fullName)

        // Check if module is already installed
        def installed = storage.getInstalledModule(reference)

        if (installed) {
            // Check integrity
            def integrity = installed.integrity
            if (integrity == ModuleIntegrity.CORRUPTED) {
                throw new AbortOperationException(
                    "Module ${reference.nameWithoutPrefix} is corrupted (missing required files). " +
                    "Please remove and reinstall."
                )
            }

            if (integrity == ModuleIntegrity.MODIFIED) {
                log.warn "Module ${reference.nameWithoutPrefix} has local modifications (checksum mismatch)"
            }

            // Check if version matches
            if (targetVersion && installed.installedVersion != targetVersion) {
                if (autoInstall) {
                    log.info "Upgrading module ${reference.nameWithoutPrefix} from ${installed.installedVersion} to ${targetVersion}"
                    return installModule(reference, targetVersion)
                } else {
                    throw new AbortOperationException(
                        "Module ${reference.nameWithoutPrefix} version mismatch: " +
                        "installed=${installed.installedVersion}, required=${targetVersion}. " +
                        "Run 'nextflow module install ${reference.nameWithoutPrefix}@${targetVersion}' to update."
                    )
                }
            }

            // Module is installed and version matches
            return installed.mainFile
        }

        // Module not installed
        if (autoInstall) {
            return installModule(reference, targetVersion)
        } else {
            throw new AbortOperationException(
                "Module ${reference.nameWithoutPrefix} is not installed. " +
                "Run 'nextflow module install ${reference.nameWithoutPrefix}' to install."
            )
        }
    }

    String resolveVersion(ModuleReference reference){
        final version = modulesConfig.getVersion(reference.fullName)
            ?: registryClient.fetchModule(reference.fullName).latest?.version
        if (!version) {
            throw new AbortOperationException("Module ${reference.nameWithoutPrefix} has no published versions")
        }
        return version
    }

    /**
     * Install or update a module
     *
     * @param reference The module reference
     * @param version Optional specific version (null = latest)
     * @param force Force reinstall even if already installed
     * @return Path to the installed module's main.nf file
     */
    Path installModule(ModuleReference reference, String version = null, boolean force = false) {
        if (!version)
            version = resolveVersion(reference)
        // Check if already installed
        if (storage.isInstalled(reference)) {
            def installed = storage.getInstalledModule(reference)
            if (installed.installedVersion == version) {
                log.info "Module ${reference.nameWithoutPrefix}@${installed.installedVersion} is already installed (version $version)"
                return installed.mainFile
            }

            // No desired version, check for local modifications
            def integrity = installed.integrity
            if (integrity == ModuleIntegrity.MODIFIED && !force) {
                throw new AbortOperationException(
                    "Module ${reference.nameWithoutPrefix} has local modifications. " +
                    "Use --force to override, or save your changes first."
                )
            }
        }


        log.info "Installing module ${reference.nameWithoutPrefix}@${version}..."

        // Download module package to temporary location
        Path tempFile = Files.createTempFile("nf-module-", ".tgz")
        try {
            // Download and validate integrity using server checksum
            def downloadResult = registryClient.downloadModule(reference.fullName, version, tempFile)

            // Install to modules directory (will compute directory checksum for future integrity checks)
            InstalledModule installed = storage.installModule(reference, version, tempFile)

            log.info "Module ${reference.nameWithoutPrefix}@${version} installed successfully at ${installed.mainFile.parent}"
            return installed.mainFile
        }
        finally {
            // Clean up temporary file
            if (Files.exists(tempFile)) {
                Files.delete(tempFile)
            }
        }
    }

}
