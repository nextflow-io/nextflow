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

import java.nio.file.Path

import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import groovy.transform.CompileStatic
import io.seqera.npr.client.RegistryClient
import nextflow.Const
import nextflow.cli.CmdRun
import nextflow.config.ConfigBuilder
import nextflow.config.RegistryConfig
import nextflow.exception.AbortOperationException
import nextflow.module.ModuleReference
import nextflow.module.ModuleResolver
import nextflow.module.RegistryClientFactory
import nextflow.util.TestOnly

/**
 * Module run subcommand
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@CompileStatic
@Parameters(commandDescription = "Run a module directly from the registry")
class CmdModuleRun extends CmdRun {
    @Parameter(names = ["-version"], description = "Module version")
    String version

    @TestOnly
    protected Path root

    @TestOnly
    protected RegistryClient client

    @Override
    String getName() {
        return 'run'
    }

    @Override
    void run() {
        final moduleFile = resolveModuleSource()
        if( moduleFile ) {
            args[0] = moduleFile.toAbsolutePath().toString()
            super.run()
        }
    }

    protected Path resolveModuleSource() {
        if( !args ) {
            throw new AbortOperationException("Module name/path not provided")
        }
        return isLocalModule(args[0])
            ? resolveLocalModule(args[0])
            : resolveRemoteModule(args[0], version)
    }

    private boolean isLocalModule(String str) {
        return str.startsWith('/') || str.startsWith('./') || str.startsWith('../')
    }

    protected Path resolveLocalModule(String str) {
        final module = Path.of(str).toAbsolutePath().normalize()
        final path = module.isDirectory() ? module.resolve(Const.DEFAULT_MAIN_FILE_NAME) : module
        if( !path.exists() )
            throw new AbortOperationException("Invalid module path: ${str}")
        return path
    }

    protected Path resolveRemoteModule(String name, String version) {
        // Parse and validate module reference
        ModuleReference reference
        try {
            reference = ModuleReference.parse(name)
        } catch( Exception e ) {
            throw new AbortOperationException("Invalid module reference: ${name}", e)
        }

        // Get config
        final baseDir = root ?: Path.of('.').toAbsolutePath().normalize()
        final config = new ConfigBuilder()
            .setOptions(launcher.options)
            .setBaseDir(baseDir)
            .build()

        final registryConfig = new RegistryConfig(config.registry as Map ?: Collections.emptyMap())
        try {
            final resolver = new ModuleResolver(baseDir, client ?: RegistryClientFactory.forConfig(registryConfig))
            return resolver.installModule(reference, version)
        } catch( Exception e ) {
            throw new AbortOperationException("Unable to install module: ${name}", e)
        }
    }
}
