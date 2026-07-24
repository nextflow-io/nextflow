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
     * Return the declared {@code requires.modules} dependencies that cannot be resolved in the
     * registry (i.e. do not exist at their pinned version). Used to validate a module's declared
     * dependencies before publishing -- the dependencies are re-resolved from the registry by
     * consumers at install time, so they must exist there (they are not bundled with the module).
     *
     * @param deps the declared dependency references ({@code scope/name@version})
     * @return the subset of {@code deps} that could not be resolved in the registry
     */
    List<String> findMissingDependencies(List<String> deps) {
        final missing = new ArrayList<String>()
        for( final dep : deps ) {
            final parsed = parseDependency(dep)
            if( !dependencyExists(parsed.reference.fullName, parsed.version) )
                missing.add(dep)
        }
        return missing
    }

    private boolean dependencyExists(String name, String version) {
        try {
            return version
                ? registryClient.getModuleRelease(name, version) != null
                : registryClient.getModule(name) != null
        }
        catch( Exception e ) {
            log.debug "Registry lookup failed for ${name}${version ? "@${version}" : ''}: ${e.message}"
            return false
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

    /**
     * @return true if the given module is installed at this resolver's base directory.
     */
    boolean isInstalled(ModuleReference reference) {
        return storage.isInstalled(reference)
    }

    /**
     * Update the vendored dependencies of an already-installed module so that its nested
     * {@code modules/} directory matches the module's (possibly locally-edited) meta.yml
     * {@code requires.modules}, without reinstalling the module itself.
     *
     * This automates the manual workflow of editing a module's meta.yml dependency versions and
     * then re-vendoring them. The parent module is left untouched -- in particular its checksum is
     * not refreshed -- so its (unpublished) modification status is preserved. Each declared
     * dependency is installed/updated at its pinned version (transitively); a locally-modified
     * dependency is not silently overwritten (an error is raised instead). A vendored dependency
     * that is no longer declared is pruned, unless it has local modifications, in which case an
     * error is raised rather than discarding the changes.
     *
     * @param reference the installed module whose dependencies should be updated
     */
    void updateDependencies(ModuleReference reference) {
        final installed = storage.getInstalledModule(reference)
        if( installed == null )
            throw new AbortOperationException("Module ${reference} is not installed")

        // declared direct dependencies from the installed (possibly locally-edited) meta.yml
        final deps = ModuleSpecFactory.fromYaml(installed.manifestFile).requiresModules
        final nestedStore = storage.nestedFor(reference)

        // remove vendored dependencies that are no longer declared
        pruneOrphanDependencies(nestedStore, deps)

        // install/update each declared dependency at its pinned version (transitively)
        final onStack = new LinkedHashSet<String>()
        onStack.add(reference.fullName)
        for( final dep : deps ) {
            final parsed = parseDependency(dep)
            walkDependencies(nestedStore, parsed.reference, parsed.version, false, onStack)
        }
    }

    /**
     * Remove vendored dependencies under {@code nestedStore} that are not in {@code declaredDeps}.
     * A modified orphan is not removed -- an error is raised so the user can reconcile it.
     */
    private static void pruneOrphanDependencies(ModuleStorage nestedStore, List<String> declaredDeps) {
        final declared = new HashSet<String>()
        for( final dep : declaredDeps )
            declared.add(parseDependency(dep).reference.fullName)

        for( final vendored : nestedStore.listInstalled() ) {
            if( vendored.reference.fullName in declared )
                continue
            if( vendored.integrity == ModuleIntegrity.MODIFIED )
                throw new AbortOperationException(
                    "Vendored dependency '${vendored.reference}' has local modifications and is no longer " +
                        "declared in meta.yml -- refusing to remove it. Restore it in meta.yml or discard the changes.")
            log.info "Removing dependency no longer declared in meta.yml: ${vendored.reference}"
            nestedStore.removeModule(vendored.reference, true)
        }
    }

    private Path walkDependencies(ModuleStorage store, ModuleReference reference, String version, boolean force,
                                  Set<String> onStack) {
        // cycle detection keys on the module name only (ignoring version): this guarantees
        // termination on the finite set of module names. A diamond (same module in sibling
        // branches) is fine — the name is popped from onStack when a branch completes.
        final key = reference.fullName

        if( key in onStack )
            throw new AbortOperationException("Module dependency cycle detected: ${(onStack.toList() << key).join(' -> ')}")

        // whether the module is already present at the requested version -- installInto will reuse
        // it as-is, so we must not touch its checksum (preserving any local modification status)
        final existing = store.getInstalledModule(reference)
        final reused = existing != null && (!version || existing.installedVersion == version)

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

        // if this module was (re)installed in this run, refresh its checksum now that its nested
        // dependencies are vendored, so the stored checksum covers the full subtree and a clean
        // install reads as VALID (the checksum is saved by installModule *before* deps are present)
        if( !reused )
            store.refreshChecksum(reference)

        return installed.mainFile
    }

}
