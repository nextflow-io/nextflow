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
import nextflow.cli.CmdRun
import nextflow.config.ConfigBuilder
import nextflow.config.ModulesConfig
import nextflow.config.RegistryConfig
import nextflow.exception.AbortOperationException
import nextflow.module.ModuleReference
import nextflow.module.ModuleResolver
import nextflow.util.NextflowSpecFile

import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths

/**
 * Module run subcommand
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Slf4j
@CompileStatic
@Parameters(commandDescription = "Run a module directly from the registry")
class ModuleRun extends CmdRun {
    @Parameter(names = ["-version"], description = "Module version")
    String version

    @Override
    String getName() {
        return 'run'
    }

    @Override
    void run() {
        if (!args ) {
            throw new AbortOperationException("Arguments not provided")
        }

        // Parse module reference (first argument starting with @)
        String moduleRef = '@' + args[0]

        // Parse and validate module reference
        ModuleReference reference
        try {
            reference = ModuleReference.parse(moduleRef)
        } catch (Exception e) {
            throw new AbortOperationException("Invalid module reference: ${moduleRef}", e)
        }

        // Get config
        def baseDir = Paths.get('.').toAbsolutePath().normalize()
        def config = new ConfigBuilder()
                .setOptions(launcher.options)
                .setBaseDir(baseDir)
                .build()

        def registryConfig = config.navigate('registry') as RegistryConfig

        //TODO: Decide final location of modules currently in nextflow_spec.json.
        // Alternative: Use nextflow config. It requires to implement nextflow.config updater features
        // def modulesConfig = config.navigate('modules') as ModulesConfig
        def specFile = new NextflowSpecFile(baseDir)
        def modulesConfig = new ModulesConfig(specFile.getModules())

        //TODO: Decide if create resolver with a temporarily storage or use current ./modules
        def tempDir = Files.createTempDirectory("nf-module-run-")
        def resolver = new ModuleResolver(tempDir, modulesConfig, registryConfig)
        try{
            Path moduleFile = resolver.installModule(reference, version)
            if( moduleFile ) {
                println "Executing module..."
                args[0] = moduleFile.toAbsolutePath().toString()
                super.run()
            }
        }
        catch (AbortOperationException e) {
            throw e
        }
        catch (Exception e) {
            log.error("Failed to run module", e)
            throw new AbortOperationException("Module run failed: ${e.message}", e)
        }
        finally {
            // Clean up temporary directory
            if (tempDir && Files.exists(tempDir)) {
                try {
                    deleteDirectory(tempDir)
                    log.debug "Cleaned up temporary directory: ${tempDir}"
                } catch (Exception e) {
                    log.warn "Failed to clean up temporary directory: ${tempDir}", e
                }
            }
        }
    }

    /**
     * Delete a directory recursively
     *
     * @param dir The directory to delete
     */
    private void deleteDirectory(Path dir) {
        if (!Files.exists(dir)) {
            return
        }

        Files.walk(dir)
                .sorted(Comparator.reverseOrder())
                .each { Path path -> Files.delete(path) }
    }
}
