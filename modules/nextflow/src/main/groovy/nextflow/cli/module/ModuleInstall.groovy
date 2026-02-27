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
import nextflow.cli.CmdBase
import nextflow.config.ConfigBuilder
import nextflow.config.ModulesConfig
import nextflow.config.RegistryConfig
import nextflow.exception.AbortOperationException
import nextflow.module.ModuleReference
import nextflow.module.ModuleRegistryClient
import nextflow.module.ModuleResolver
import nextflow.pipeline.PipelineSpec
import nextflow.util.TestOnly

import java.nio.file.Path
import java.nio.file.Paths

/**
 * Module install subcommand
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Slf4j
@Parameters(commandDescription = "Install a module from the registry")
@CompileStatic
class ModuleInstall extends CmdBase {

    @Parameter(names = ["-version"], description = "Module version")
    String version

    @Parameter(names = ["-force"], description = "Force reinstall even if already installed", arity = 0)
    boolean force = false

    @Parameter(description = "[scope/name]", required = true)
    List<String> args

    @TestOnly
    protected Path root

    @TestOnly
    protected ModuleRegistryClient client

    @Override
    String getName() {
        return 'install'
    }

    @Override
    void run() {
        if (!args || args.size() != 1) {
            throw new AbortOperationException("Incorrect number of arguments")
        }

        def moduleRef = '@' + args[0]

        def reference = ModuleReference.parse(moduleRef)

        // Get config
        def baseDir = root ?: Paths.get('.').toAbsolutePath().normalize()
        def config = new ConfigBuilder()
                .setOptions(launcher.options)
                .setBaseDir(baseDir)
                .build()
        final registryConfig = config.navigate('registry') as RegistryConfig ?: new RegistryConfig()

        // Get modules versions from nextflow_spec.json.
        final specFile = new PipelineSpec(baseDir)
        final modulesConfig = new ModulesConfig(specFile.getModules())

        // Create resolver and install
        def resolver = new ModuleResolver(baseDir, client ?: new ModuleRegistryClient(registryConfig), modulesConfig)

        try {
            def installedMainFile = resolver.installModule(reference, version, force)

            // Update nextflow_spec.json with the installed module version
            def installedVersion = version ?: resolver.resolveVersion(reference)
            specFile.addModuleEntry(reference.fullName, installedVersion)

            println "Module ${reference.nameWithoutPrefix}@${installedVersion} installed and configured successfully"
        }
        catch (AbortOperationException e) {
            throw e
        }
        catch (Exception e) {
            throw new AbortOperationException("Installation failed: ${e.message}", e)
        }
    }
}
