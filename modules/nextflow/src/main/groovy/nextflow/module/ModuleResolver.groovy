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
import io.seqera.npr.client.RegistryClient

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

    private final RegistryClient registryClient
    private final ModuleStorage storage

    ModuleResolver(Path baseDir, RegistryClient registryClient) {
        this.registryClient = registryClient
        this.storage = new ModuleStorage(baseDir)
    }

    ModuleResolver(Path baseDir, RegistryConfig registryConfig = null) {
        this(baseDir, RegistryClientFactory.forConfig(registryConfig ?: new RegistryConfig()))

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

        // Check if module is already installed
        def installed = storage.getInstalledModule(reference)

        if( installed ) {
            // Check integrity
            def integrity = installed.integrity
            if( integrity == ModuleIntegrity.CORRUPTED ) {
                throw new AbortOperationException(
                    "Module ${reference} is corrupted (missing required files). " +
                        "Please remove and reinstall."
                )
            }

            if( integrity == ModuleIntegrity.MODIFIED ) {
                log.warn1 "Module ${reference} has local modifications (checksum mismatch)"
            } else if( integrity == ModuleIntegrity.NO_REMOTE_MODULE ) {
                log.warn1 "Module ${reference} has no registry origin (.module-info missing)"
            }

            // Check if version matches
            if( version && installed.installedVersion != version ) {
                if( autoInstall ) {
                    log.info "Upgrading module ${reference} from ${installed.installedVersion} to ${version}"
                    return installWithDependencies(reference, version)
                } else {
                    throw new AbortOperationException(
                        "Module ${reference} version mismatch: " +
                            "installed=${installed.installedVersion}, required=${version}. " +
                            "Run 'nextflow module install ${reference}@${version}' to update."
                    )
                }
            }

            // Module is installed and version matches
            return installed.mainFile
        }

        // Module not installed
        if( autoInstall ) {
            return installWithDependencies(reference, version)
        } else {
            throw new AbortOperationException(
                "Module ${reference} is not installed. " +
                    "Run 'nextflow module install ${reference}' to install."
            )
        }
    }

    String resolveVersion(ModuleReference reference) {
        final version = registryClient.getModule(reference.fullName)?.latest?.version
        if( !version ) {
            throw new AbortOperationException("Module ${reference} has no published versions")
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
        return installInto(storage, reference, version, force).mainFile
    }

    /**
     * Install or update a module into the given storage and return the resulting
     * {@link InstalledModule}. The storage's base directory determines where the
     * module is placed -- the project root for a top-level module, or a parent
     * module's own directory for a nested (vendored) dependency.
     *
     * @param store     the target storage (defines the install base directory)
     * @param reference the module reference
     * @param version   specific version (null = latest)
     * @param force     force reinstall even if locally modified
     * @return the installed module
     */
    private InstalledModule installInto(ModuleStorage store, ModuleReference reference, String version, boolean force) {
        // Check if already installed locally before hitting the registry
        if( store.isInstalled(reference) ) {
            def installed = store.getInstalledModule(reference)

            // No specific version requested -- use the local module as-is
            if( !version ) {
                log.debug "Module ${reference}@${installed.installedVersion} is already installed locally"
                return installed
            }

            if( installed.installedVersion == version ) {
                log.debug "Module ${reference}@${installed.installedVersion} is already installed (version $version)"
                return installed
            }

            // Version mismatch -- check for local modifications before overwriting
            def integrity = installed.integrity
            if( integrity == ModuleIntegrity.MODIFIED && !force ) {
                throw new AbortOperationException(
                    "Module ${reference} has local modifications. " +
                        "Use '-force' to override, or save your changes first."
                )
            }
            if( integrity == ModuleIntegrity.NO_REMOTE_MODULE && !force ) {
                throw new AbortOperationException(
                    "Folder 'modules/${reference}' already exists and is not a valid remote module. " +
                        "Use '-force' to override, or save your changes first."
                )
            }
        }

        if( !version )
            version = resolveVersion(reference)

        log.info "Installing module ${reference}@${version}..."

        // Download module package to temporary location
        Path tempFile = Files.createTempFile("nf-module-", ".tgz")
        try {
            // Download and validate integrity using server checksum
            def downloadUrl = registryClient.downloadModuleRelease(reference.fullName, version, tempFile)

            // Install to modules directory (will compute directory checksum for future integrity checks)
            InstalledModule installed = store.installModule(reference, version, tempFile, downloadUrl)

            log.info "Module ${reference}@${version} installed successfully at ${installed.mainFile.parent.toAbsolutePath()}"
            return installed
        }
        finally {
            // Clean up temporary file
            if( Files.exists(tempFile) ) {
                Files.delete(tempFile)
            }
        }
    }

    /**
     * A parsed {@code requires.modules} entry: a module reference plus an
     * optional pinned version (the part after {@code @}).
     */
    @groovy.transform.TupleConstructor
    static class DependencySpec {
        ModuleReference reference
        String version
    }

    /**
     * Parse a {@code requires.modules} entry of the form
     * {@code scope/name[@<version-or-constraint>]}.
     */
    protected static DependencySpec parseDependency(String dep) {
        final idx = dep.lastIndexOf('@')
        if( idx < 0 )
            return new DependencySpec(ModuleReference.parse(dep), null)
        final name = dep.substring(0, idx)
        final version = dep.substring(idx + 1)
        return new DependencySpec(ModuleReference.parse(name), version ?: null)
    }

    /**
     * Install a module together with its transitive {@code requires.modules}
     * dependencies, using nested per-module vendoring (ADR v2, PR #7342).
     *
     * The module installs under this resolver's base (`<base>/modules/<scope>/<name>/`);
     * each of its *direct* dependencies is vendored under the module's own nested
     * `modules/` directory, and each dependency likewise vendors its own
     * dependencies (arbitrary depth). Dependencies are installed at their pinned
     * version in isolation -- there is no cross-module flattening or conflict
     * resolution, so the same module may be vendored more than once in the tree.
     * Dependency cycles are detected and reported rather than followed.
     *
     * @param reference the root module reference
     * @param version   optional explicit version for the root (null = latest)
     * @param force     force reinstall of locally modified modules
     * @return the path to the root module's main.nf file
     */
    Path installWithDependencies(ModuleReference reference, String version = null, boolean force = false) {
        return walkDependencies(storage, reference, version, force, new LinkedHashSet<String>())
    }

    private Path walkDependencies(ModuleStorage store, ModuleReference reference, String version, boolean force,
                                  Set<String> onStack) {
        // cycle detection keys on the module name only (ignoring version): this guarantees
        // termination on the finite set of module names. A diamond (same module in sibling
        // branches) is fine — the name is popped from onStack when a branch completes.
        final key = reference.fullName

        if( key in onStack )
            throw new AbortOperationException("Module dependency cycle detected: ${(onStack.toList() << key).join(' -> ')}")

        // install this module at the current level
        final installed = installInto(store, reference, version, force)

        // read its direct dependencies from the installed spec (requiresModules is never null)
        final deps = ModuleSpecFactory.fromYaml(installed.manifestFile).requiresModules

        // dependencies are vendored under THIS module's own nested `modules/` directory
        final nestedStore = store.nestedFor(reference)
        onStack.add(key)
        for( final dep : deps ) {
            final parsed = parseDependency(dep)
            walkDependencies(nestedStore, parsed.reference, parsed.version, force, onStack)
        }
        onStack.remove(key)

        return installed.mainFile
    }

}
