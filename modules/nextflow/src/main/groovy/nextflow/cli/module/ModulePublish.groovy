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

package nextflow.cli.module

import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Const
import nextflow.cli.CmdBase
import nextflow.config.ConfigBuilder
import nextflow.config.RegistryConfig
import nextflow.exception.AbortOperationException
import nextflow.module.ModuleManifest
import nextflow.module.ModuleReference
import nextflow.module.ModuleRegistryClient
import nextflow.module.ModuleStorage

import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths

/**
 * Module publish subcommand
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Slf4j
@CompileStatic
@Parameters(commandDescription = "Publish a module to the registry")
class ModulePublish extends CmdBase {

    @Parameter(names = ["-dry-run"], description = "Validate without uploading", arity=0)
    boolean dryRun = false

    @Parameter(names = ["-registry"], description = "Target registry URL")
    String registryUrl = RegistryConfig.DEFAULT_REGISTRY_URL

    @Parameter(description = "Module directory path or scope/name")
    List<String> args

    @Override
    String getName() {
        return 'publish'
    }

    @Override
    void run() {
        if (!args || args.size() != 1) {
            throw new AbortOperationException("Incorrect number of arguments")
        }

        Path moduleDir = determineModuleDir(args[0])

        log.info "Publishing module from: ${moduleDir}"

        // Step 1: Validate module structure
        def validationErrors = validateModuleStructure(moduleDir)
        if (!validationErrors.isEmpty()) {
            throw new AbortOperationException(
                "Module validation failed:\n" + validationErrors.collect { "  - ${it}" }.join('\n')
            )
        }

        // Step 2: Load and validate manifest
        def manifestPath = moduleDir.resolve(ModuleStorage.MODULE_MANIFEST_FILE)
        def manifest = ModuleManifest.load(manifestPath)

        def manifestErrors = manifest.validate()
        if (!manifestErrors.isEmpty()) {
            throw new AbortOperationException(
                "Module manifest validation failed:\n" + manifestErrors.collect { "  - ${it}" }.join('\n')
            )
        }

        log.info "Module validated: ${manifest.name}@${manifest.version}"

        if (dryRun) {
            printDryRunInfo(manifest)
            return
        }

        // Step 3: Get authentication token
        def config = new ConfigBuilder()
            .setOptions(launcher.options)
            .setBaseDir(moduleDir)
            .build()

        def registryConfig = config.navigate('registry') as RegistryConfig

        publishModule(moduleDir, registryConfig, manifest)

    }

    private void publishModule(Path moduleDir, RegistryConfig registryConfig, ModuleManifest manifest){
        log.info "Creating module bundle..."
        def storage = new ModuleStorage(moduleDir.parent)
        def tempBundleFile = Files.createTempFile("nf-module-publish-", ".tar.gz")

        try {
            storage.createBundle(moduleDir, tempBundleFile)

            // Compute bundle checksum
            def checksum = storage.computeBundleChecksum(tempBundleFile)
            log.info "Bundle checksum: ${checksum}"

            // Read bundle content as bytes
            def bundleBytes = Files.readAllBytes(tempBundleFile)

            // Create publish request as a map (npr-api will serialize it)
            def request = [
                version: manifest.version,
                bundle: bundleBytes
            ]

            // Publish to registry
            log.info "Publishing module to registry: ${registryUrl}"
            def registryClient = new ModuleRegistryClient(registryConfig)
            def response = registryClient.publishModule(manifest.name, request, registryUrl)

            println "✓ Module published successfully!"
            println ""
            println "Module details:"
            println "  Name: ${manifest.name}"
            println "  Version: ${manifest.version}"
            println "  DownloadUrl: ${response.downloadUrl}"

            println ""
            println "Others can now install this module using:"
            println "  nextflow module install ${manifest.name}"

        } finally {
            // Clean up temporary bundle file
            if (Files.exists(tempBundleFile)) {
                try {
                    Files.delete(tempBundleFile)
                } catch (Exception e) {
                    log.warn "Failed to clean up temporary bundle file: ${e.message}"
                }
            }
        }
    }

    private void printDryRunInfo(ModuleManifest manifest) {
        println "✓ Module structure is valid"
        println ""
        println "Module details:"
        println "  Name: ${manifest.name}"
        println "  Version: ${manifest.version}"
        println "  Description: ${manifest.description}"
        println "  License: ${manifest.license}"
        if( manifest.authors ) {
            println "  Authors: ${manifest.authors.join(', ')}"
        }
        if( manifest.keywords ) {
            println "  Keywords: ${manifest.keywords.join(', ')}"
        }
        if( manifest.requires ) {
            println "  Requires:"
            manifest.requires.each { name, version ->
                println "    - ${name}: ${version}"
            }
        }
        println ""
        println "Dry run complete. Module is ready to publish."
        println "Run without --dry-run to publish to the registry."
    }

    /**
     * Validate that the module directory has the required structure
     *
     * @param moduleDir The module directory path
     * @return List of validation error messages (empty if valid)
     */
    private List<String> validateModuleStructure(Path moduleDir) {
        List<String> errors = []

        if (!Files.exists(moduleDir) || !Files.isDirectory(moduleDir)) {
            errors << "Module directory does not exist: ${moduleDir}".toString()
            return errors
        }

        // Check for required files
        def mainNf = moduleDir.resolve(Const.DEFAULT_MAIN_FILE_NAME)
        if (!Files.exists(mainNf)) {
            errors << "Missing required file: $Const.DEFAULT_MAIN_FILE_NAME".toString()
        }

        def metaYaml = moduleDir.resolve(ModuleStorage.MODULE_MANIFEST_FILE)
        if (!Files.exists(metaYaml)) {
            errors << "Missing required file: $ModuleStorage.MODULE_MANIFEST_FILE".toString()
        }

        def readme = moduleDir.resolve(ModuleStorage.MODULE_README_FILE)
        if (!Files.exists(readme)) {
            errors << "Missing required file: $ModuleStorage.MODULE_README_FILE".toString()
        }

        // Check bundle size (1MB uncompressed limit)
        try {
            long totalSize = Files.walk(moduleDir)
                .filter { Files.isRegularFile(it) }
                .mapToLong { Files.size(it) }
                .sum()

            def maxSize = 1024 * 1024 // 1MB in bytes
            if (totalSize > maxSize) {
                def sizeMB = totalSize / (1024 * 1024)
                errors << "Module size exceeds 1MB limit (current: ${String.format('%.2f', sizeMB)}MB)".toString()
            }
        } catch (Exception e) {
            log.warn "Failed to check module size: ${e.message}"
        }

        return errors
    }
    /**
     * Determine if the specified module is a local path or a reference
     * @param module
     * @return
     */
    private Path determineModuleDir(String module) {
        //If local path exists return this path as module dir
        if (Paths.get(module).exists()){
            return Paths.get(module).toAbsolutePath().normalize()
        }

        final ref =  ModuleReference.parse('@' + module)
        final localStorage = new ModuleStorage(Paths.get('.').toAbsolutePath().normalize())

        if (!localStorage.isInstalled(ref)){
            throw new AbortOperationException("No module diretory found for $module")
        }

        return localStorage.getModuleDir(ref)
    }
}
