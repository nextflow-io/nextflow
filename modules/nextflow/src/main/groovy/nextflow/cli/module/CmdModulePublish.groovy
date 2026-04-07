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
import io.seqera.npr.client.RegistryClient
import nextflow.cli.CmdBase
import nextflow.config.ConfigBuilder
import nextflow.config.RegistryConfig
import nextflow.exception.AbortOperationException
import nextflow.module.ModuleChecksum
import nextflow.module.ModuleInfo
import nextflow.module.ModuleSpec
import nextflow.module.ModuleSpecFactory
import nextflow.module.ModuleReference
import nextflow.module.ModuleValidator
import nextflow.module.RegistryClientFactory
import nextflow.module.ModuleStorage
import nextflow.util.TestOnly

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
class CmdModulePublish extends CmdBase {

    @Parameter(names = ["-dry-run"], description = "Validate without uploading", arity=0)
    boolean dryRun = false

    @Parameter(names = ["-registry"], description = "Target registry URL.")
    String registryUrl

    @Parameter(description = "Module directory path or scope/name")
    List<String> args

    @TestOnly
    protected Path root

    @TestOnly
    protected RegistryClient client

    //Flag if publish is invoked from a scope/name. In this case we should create/update the .module-info with the correct checksum
    private boolean useModuleReference = false

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

        // Step 1: Validate module structure and spec
        def validationErrors = ModuleValidator.validate(moduleDir)
        if (!validationErrors.isEmpty()) {
            throw new AbortOperationException(
                "Module validation failed:\n" + validationErrors.collect { "  - ${it}" }.join('\n')
            )
        }

        // Step 2: Load spec for publish metadata
        def manifestPath = moduleDir.resolve(ModuleStorage.MODULE_MANIFEST_FILE)
        def spec = ModuleSpecFactory.fromYaml(manifestPath)

        log.info "Module validated: ${spec.name}@${spec.version}"

        if (dryRun) {
            printDryRunInfo(spec)
            return
        }

        // Step 3: Get authentication token
        def config = new ConfigBuilder()
            .setOptions(launcher.options)
            .setBaseDir(moduleDir)
            .build()

        def registryConfig = config.navigate('registry') as RegistryConfig ?: new RegistryConfig()

        publishModule(moduleDir, registryConfig, spec)

    }

    private void publishModule(Path moduleDir, RegistryConfig registryConfig, ModuleSpec spec){
        log.info "Creating module bundle..."
        def tempBundleFile = Files.createTempFile("nf-module-publish-", ".tar.gz")

        try {
            def checksum = ModuleStorage.createBundle(moduleDir, tempBundleFile)
            log.info "Bundle checksum: ${checksum}"

            // Read bundle content as bytes
            def bundleBytes = Files.readAllBytes(tempBundleFile)

            // Create publish request as a map (npr-api will serialize it)
            def request = [
                version: spec.version,
                bundle: bundleBytes
            ]

            // Publish to registry
            final registry = registryUrl ?: registryConfig.url
            log.info "Publishing module to registry: ${registryUrl ?: registryConfig.url}"
            def registryClient = RegistryClientFactory.forConfig(registryConfig)
            def response = registryClient.publishModuleRelease(spec.name, request, registry)

            if (useModuleReference) {
                // If publish is performed using the module reference we should create/update the .module-info with the correct checksum
                try {
                    ModuleInfo.save(moduleDir, [checksum: ModuleChecksum.compute(moduleDir), registryUrl: registry] )
                }catch (Exception e){
                    log.warn("Unable to save the checksum - ${e.message}")
                }
            }
            println "✓ Module published successfully!"
            println ""
            println "Module details:"
            println "  Name: ${spec.name}"
            println "  Version: ${spec.version}"
            println "  DownloadUrl: ${response.downloadUrl}"

            println ""
            println "Others can now install this module using:"
            println "  nextflow module install ${spec.name}"

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

    private void printDryRunInfo(ModuleSpec spec) {
        println "✓ Module structure is valid"
        println ""
        println "Module details:"
        println "  Name: ${spec.name}"
        println "  Version: ${spec.version}"
        println "  Description: ${spec.description}"
        println "  License: ${spec.license}"
        if( spec.authors ) {
            println "  Authors: ${spec.authors.join(', ')}"
        }
        if( spec.keywords ) {
            println "  Keywords: ${spec.keywords.join(', ')}"
        }
        if( spec.requires ) {
            println "  Requires:"
            spec.requires.each { name, version ->
                println "    - ${name}: ${version}"
            }
        }
        println ""
        println "Dry run complete. Module is ready to publish."
        println "Run without -dry-run to publish to the registry."
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

        final ref =  ModuleReference.parse(module)
        final localStorage = new ModuleStorage(root ?: Paths.get('.').toAbsolutePath().normalize())

        if (!localStorage.isInstalled(ref)){
            throw new AbortOperationException("No module directory found for $module")
        }
        useModuleReference = true
        return localStorage.getModuleDir(ref)
    }
}
