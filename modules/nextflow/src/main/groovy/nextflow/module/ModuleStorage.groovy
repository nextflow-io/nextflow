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
import nextflow.Const
import nextflow.exception.AbortOperationException
import nextflow.file.FileHelper
import org.apache.commons.compress.archivers.tar.TarArchiveEntry
import org.apache.commons.compress.archivers.tar.TarArchiveInputStream
import org.apache.commons.compress.archivers.tar.TarArchiveOutputStream
import org.apache.commons.compress.compressors.gzip.GzipCompressorInputStream
import org.apache.commons.compress.compressors.gzip.GzipCompressorOutputStream

import java.nio.file.Files
import java.nio.file.Path
import java.util.stream.Stream
import java.util.zip.ZipEntry
import java.util.zip.ZipInputStream

import static nextflow.module.ModuleInfo.MODULE_INFO_FILE

/**
 * Manages local filesystem storage for modules
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Slf4j
@CompileStatic
class ModuleStorage {
    public static final String MODULE_MANIFEST_FILE = "meta.yml"
    public static final String MODULE_README_FILE = "README.md"
    private final Path modulesDir

    /**
     * Create a ModuleStorage instance
     *
     * @param baseDir The base directory (usually project root)
     */
    ModuleStorage(Path baseDir) {
        this.modulesDir = baseDir.resolve('modules')
    }

    /**
     * Get the modules directory path
     *
     * @return The modules directory
     */
    Path getModulesDir() {
        return modulesDir
    }

    /**
     * Get the directory path for a specific module
     *
     * @param reference The module reference
     * @return The module directory path
     */
    Path getModuleDir(ModuleReference reference) {
        return modulesDir.resolve(reference.scope).resolve(reference.name)
    }

    /**
     * Get the module info path for a specific module
     *
     * @param reference The module reference
     * @return The module info path
     */
    Path getModuleInfo(ModuleReference reference) {
        return modulesDir.resolve(reference.scope).resolve(reference.name).resolve(MODULE_INFO_FILE)
    }

    /**
     * Check if a module is installed locally
     *
     * @param reference The module reference
     * @return true if the module directory exists
     */
    boolean isInstalled(ModuleReference reference) {
        def moduleDir = getModuleDir(reference)
        return Files.exists(moduleDir) && Files.isDirectory(moduleDir)
    }

    /**
     * Get an installed module
     *
     * @param reference The module reference
     * @return InstalledModule object, or null if not installed
     */
    InstalledModule getInstalledModule(ModuleReference reference) {
        def moduleDir = getModuleDir(reference)
        if (!Files.exists(moduleDir) || !Files.isDirectory(moduleDir)) {
            return null
        }

        def installed = new InstalledModule(
            reference: reference,
            directory: moduleDir,
            mainFile: moduleDir.resolve(Const.DEFAULT_MAIN_FILE_NAME),
            manifestFile: moduleDir.resolve(MODULE_MANIFEST_FILE),
            moduleInfoFile: moduleDir.resolve(MODULE_INFO_FILE),
        )

        // Load checksum if available
        Map infoProps = ModuleInfo.load(moduleDir)
        installed.expectedChecksum = infoProps?.checksum
        installed.registryUrl = infoProps?.registryUrl
        installed.installedVersion = ModuleSpec.load(installed.manifestFile).version
        return installed
    }

    /**
     * List all installed modules by scanning for directories containing a .module-info marker file.
     *
     * @return List of InstalledModule objects
     */
    List<InstalledModule> listInstalled() {
        if (!Files.exists(modulesDir) || !Files.isDirectory(modulesDir)) {
            return []
        }

        List<InstalledModule> modules = []

        try( final walkStream = Files.walk(modulesDir) ) {
            walkStream
                .filter { Path path -> Files.isDirectory(path) }
                .filter { Path path -> Files.exists(path.resolve(MODULE_INFO_FILE)) }
                .each { Path moduleDir ->
                    try {
                        def rel = modulesDir.relativize(moduleDir)
                        if( rel.nameCount < 2 ) return  // Need at least scope/name
                        def reference = ModuleReference.parse(rel.toString())
                        def installed = getInstalledModule(reference)
                        if( installed ) modules.add(installed)
                    } catch(Exception e){
                        // Catching exception to go on inspecting other valid folders
                        log.debug("Not a valid module reference - $e.message")
                    }
                }
        } catch (IOException e) {
            log.warn "Failed to scan modules directory ${modulesDir}: ${e.message}"
        }

        return modules
    }

    /**
     * Install a module from a downloaded package file
     *
     * @param reference The module reference
     * @param version The module version
     * @param packageFile The downloaded package file (zip or tar.gz)
     * @return The InstalledModule object
     */
    InstalledModule installModule(ModuleReference reference, String version, Path packageFile, String downloadUrl) {
        def moduleDir = getModuleDir(reference)

        try {
            // Remove existing installation if present
            if (Files.exists(moduleDir)) {
                log.debug "Removing existing module installation: ${moduleDir}"
                FileHelper.deletePath(moduleDir)
            }

            // Create module directory
            Files.createDirectories(moduleDir)

            // Extract package - detect format by file extension
            if (packageFile.toString().endsWith('.tgz') || packageFile.toString().endsWith('.tar.gz')) {
                extractTarGz(packageFile, moduleDir)
            } else {
                extractZip(packageFile, moduleDir)
            }

            // Compute and save checksum of extracted directory contents
            // This checksum is used to detect local modifications
            def checksum = ModuleChecksum.compute(moduleDir)
            ModuleInfo.save(moduleDir, [checksum: checksum, registryUrl: downloadUrl])

            log.debug "Installed module ${reference}@${version} to ${moduleDir}"

            return getInstalledModule(reference)
        }
        catch (Exception e) {
            // Clean up on failure
            if (Files.exists(moduleDir)) {
                try {
                    FileHelper.deletePath(moduleDir)
                } catch (Exception cleanupError) {
                    log.warn "Failed to clean up after installation failure: ${cleanupError.message}"
                }
            }
            throw new AbortOperationException("Failed to install module ${reference}@${version}", e)
        }
    }

    /**
     * Remove an installed module
     *
     * @param reference The module reference
     * @param force Force local module folder
     * @return true if module was removed, false if not installed
     */
    boolean removeModule(ModuleReference reference, boolean force) {
        final installed = getInstalledModule(reference)
        if( !installed )
            return
        final integrity = installed.integrity
        if( integrity == ModuleIntegrity.NO_REMOTE_MODULE && !force ) {
            throw new AbortOperationException(
                "Folder 'modules/${reference}' already exists and is not a valid remote module ($MODULE_INFO_FILE missing). " +
                    "Use '-force' to remove, or save your changes first.")
        }

        if( integrity == ModuleIntegrity.MODIFIED && !force ) {
            throw new AbortOperationException(
                "Module ${reference} has local modifications. " +
                    "Use '-force' to remove, or save your changes first."
            )
        }

        try {
            FileHelper.deletePath(installed.directory)
            log.debug "Removed module: ${reference}"

            // Clean up empty scope directory
            def scopeDir = installed.directory.parent
            if (Files.exists(scopeDir) && isEmpty(scopeDir)) {
                Files.delete(scopeDir)
            }

            return true
        }
        catch (Exception e) {
            throw new AbortOperationException("Failed to remove module ${reference}", e)
        }
    }

    /**
     * Extract a zip file to a target directory
     *
     * @param zipFile The zip file path
     * @param targetDir The target directory
     */
    private void extractZip(Path zipFile, Path targetDir) {
        Files.newInputStream(zipFile).withCloseable { fis ->
            new ZipInputStream(fis).withCloseable { zis ->
                ZipEntry entry
                while ((entry = zis.nextEntry) != null) {
                    def targetPath = targetDir.resolve(entry.name)

                    // Security check: prevent zip slip
                    if (!targetPath.normalize().startsWith(targetDir.normalize())) {
                        throw new AbortOperationException("Invalid zip entry: ${entry.name}")
                    }

                    if (entry.directory) {
                        Files.createDirectories(targetPath)
                    } else {
                        // Create parent directories if needed
                        if (targetPath.parent) {
                            Files.createDirectories(targetPath.parent)
                        }

                        // Write file
                        Files.copy(zis, targetPath)
                    }

                    zis.closeEntry()
                }
            }
        }
    }

    /**
     * Extract a tar.gz file to a target directory
     *
     * @param tarGzFile The tar.gz file path
     * @param targetDir The target directory
     */
    private void extractTarGz(Path tarGzFile, Path targetDir) {
        Files.newInputStream(tarGzFile).withCloseable { fis ->
            new GzipCompressorInputStream(fis).withCloseable { gzis ->
                new TarArchiveInputStream(gzis).withCloseable { tis ->
                    TarArchiveEntry entry
                    while ((entry = tis.nextTarEntry) != null) {
                        def targetPath = targetDir.resolve(entry.name)

                        // Security check: prevent tar slip
                        if (!targetPath.normalize().startsWith(targetDir.normalize())) {
                            throw new AbortOperationException("Invalid tar entry: ${entry.name}")
                        }

                        if (entry.directory) {
                            Files.createDirectories(targetPath)
                        } else {
                            // Create parent directories if needed
                            if (targetPath.parent) {
                                Files.createDirectories(targetPath.parent)
                            }

                            // Write file
                            Files.copy(tis, targetPath)
                        }
                    }
                }
            }
        }
    }

    /**
     * Check if a directory is empty
     *
     * @param dir The directory to check
     * @return true if directory is empty
     */
    private boolean isEmpty(Path dir) {
        if (!Files.exists(dir) || !Files.isDirectory(dir)) {
            return true
        }

        try {
            Stream<Path> entries = Files.list(dir)
            try {
                return !entries.findFirst().isPresent()
            } finally {
                entries.close()
            }
        } catch (Exception e) {
            log.warn "Failed to check if directory is empty: ${dir}", e
            return false
        }
    }

    /**
     * Create a module bundle (tar.gz) from a module directory for publishing
     *
     * @param moduleDir The module directory to bundle
     * @param targetFile The target bundle file path
     * @return The created bundle file checksum
     */
    static String createBundle(Path moduleDir, Path targetFile) {
        if (!Files.exists(moduleDir) || !Files.isDirectory(moduleDir)) {
            throw new AbortOperationException("Module directory not found: ${moduleDir}")
        }

        try {
            // Create parent directories if needed
            if (targetFile.parent) {
                Files.createDirectories(targetFile.parent)
            }

            // Create tar.gz bundle
            Files.newOutputStream(targetFile).withCloseable { fos ->
                new GzipCompressorOutputStream(fos).withCloseable { gzos ->
                    new TarArchiveOutputStream(gzos).withCloseable { tos ->
                        // Add all files in module directory to bundle
                        addToTarArchive(tos, moduleDir, moduleDir)
                    }
                }
            }
            final checksum = computeBundleChecksum(targetFile)
            log.debug "Created module bundle: ${targetFile} (size: ${Files.size(targetFile)} bytes, checksum: $checksum)"
            return checksum
        }
        catch (Exception e) {
            // Clean up partial file on failure
            if (Files.exists(targetFile)) {
                try {
                    Files.delete(targetFile)
                } catch (Exception cleanupError) {
                    log.warn "Failed to clean up after bundle creation failure: ${cleanupError.message}"
                }
            }
            throw new AbortOperationException("Failed to create module bundle", e)
        }
    }

    /**
     * Add files to tar archive recursively
     *
     * @param tos The tar archive output stream
     * @param sourceDir The source directory being archived
     * @param currentPath The current path being added
     */
    private static void addToTarArchive(TarArchiveOutputStream tos, Path sourceDir, Path currentPath) {

        try ( def tarStream = Files.list(currentPath)) {
            tarStream.each { Path path ->
                // Skip .module-info file when creating bundle
                if (path.fileName.toString() == MODULE_INFO_FILE) {
                    return
                }

                def relativePath = sourceDir.relativize(path).toString()

                if (Files.isDirectory(path)) {
                    // Add directory entry
                    def entry = new TarArchiveEntry(path.toFile(), "${relativePath}/")
                    tos.putArchiveEntry(entry)
                    tos.closeArchiveEntry()

                    // Recursively add directory contents
                    addToTarArchive(tos, sourceDir, path)
                } else {
                    // Add file entry
                    def entry = new TarArchiveEntry(path.toFile(), relativePath)
                    entry.setSize(Files.size(path))
                    tos.putArchiveEntry(entry)

                    // Copy file content
                    Files.copy(path, tos)
                    tos.closeArchiveEntry()
                }
            }
        }
    }

    /**
     * Compute the checksum of a module bundle
     *
     * @param bundleFile The bundle file
     * @return The SHA-256 checksum as hex string
     */
    static String computeBundleChecksum(Path bundleFile) {
        if (!Files.exists(bundleFile)) {
            throw new AbortOperationException("Bundle file not found: ${bundleFile}")
        }

        try {
            return ModuleChecksum.computeFile(bundleFile)
        }
        catch (Exception e) {
            throw new AbortOperationException("Failed to compute bundle checksum", e)
        }
    }
}
