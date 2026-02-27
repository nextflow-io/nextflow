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
import nextflow.cli.CmdRun
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
 * Module run subcommand
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@CompileStatic
@Parameters(commandDescription = "Run a module directly from the registry")
class ModuleRun extends CmdRun {
    @Parameter(names = ["-version"], description = "Module version")
    String version

    @TestOnly
    protected Path root

    @TestOnly
    protected ModuleRegistryClient client

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
        def baseDir = root ?: Paths.get('.').toAbsolutePath().normalize()
        def config = new ConfigBuilder()
                .setOptions(launcher.options)
                .setBaseDir(baseDir)
                .build()

        def registryConfig = config.navigate('registry') as RegistryConfig ?: new RegistryConfig()

        //Get module version from nextflow_spec.json.
        def specFile = new PipelineSpec(baseDir)
        def modulesConfig = new ModulesConfig(specFile.getModules())

        def resolver = new ModuleResolver(baseDir, client ?: new ModuleRegistryClient(registryConfig), modulesConfig)
        Path moduleFile = resolver.installModule(reference, version)
        if( moduleFile ) {
            println "Executing module..."
            args[0] = moduleFile.toAbsolutePath().toString()
            super.run()
        }
    }
}
