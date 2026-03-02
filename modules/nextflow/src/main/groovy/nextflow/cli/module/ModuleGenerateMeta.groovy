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
import nextflow.Session
import nextflow.cli.CmdBase
import nextflow.config.ConfigBuilder
import nextflow.exception.AbortOperationException
import nextflow.module.MetaYmlGenerator
import nextflow.util.TestOnly

import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths

/**
 * Implements the {@code nextflow module generate-meta} subcommand.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Slf4j
@CompileStatic
@Parameters(commandDescription = "Generate a meta.yml manifest for a local module")
class ModuleGenerateMeta extends CmdBase {

    @Parameter(description = "Module directory path")
    List<String> args

    @Parameter(names = ["-force"], description = "Overwrite existing meta.yml", arity = 0)
    boolean force = false

    @Parameter(names = ["-dry-run"], description = "Print generated YAML to stdout without writing any file", arity = 0)
    boolean dryRun = false

    @Parameter(names = ["-name"], description = "Module name in 'scope/name' format (e.g. nf-core/fastqc)")
    String moduleName

    @Parameter(names = ["-version"], description = "Module version (e.g. 1.0.0)")
    String moduleVersion

    @Parameter(names = ["-description"], description = "Short description of what the module does")
    String description

    @Parameter(names = ["-license"], description = "Module license identifier (e.g. MIT, Apache-2.0)")
    String license

    @Parameter(names = ["-author"], description = "Module author; may be specified multiple times")
    List<String> authors

    @TestOnly
    protected Path root

    @Override
    String getName() { 'generate-meta' }

    @Override
    void run() {
        final baseDir = root ?: Paths.get('.').toAbsolutePath().normalize()
        final moduleDir = resolveModuleDir(baseDir)
        final mainNfPath = moduleDir.resolve('main.nf')

        if( !Files.exists(mainNfPath) )
            throw new AbortOperationException("Missing required file: main.nf in ${moduleDir}")

        final metadata = parseModule(mainNfPath)
        final renderOptions = new MetaYmlGenerator.RenderOptions(
            name: moduleName,
            version: moduleVersion,
            description: description,
            license: license,
            authors: authors ?: []
        )
        final yamlContent = MetaYmlGenerator.render(metadata, renderOptions)

        if( dryRun ) {
            println yamlContent
            return
        }

        final metaYmlPath = moduleDir.resolve('meta.yml')
        if( Files.exists(metaYmlPath) && !force )
            throw new AbortOperationException("meta.yml already exists. Use -force to overwrite.")

        Files.writeString(metaYmlPath, yamlContent)
        println "Generated: ${metaYmlPath}"
    }

    private Path resolveModuleDir(Path baseDir) {
        if( !args ) {
            return baseDir
        }
        final path = Paths.get(args[0])
        final moduleDir = path.absolute ? path : baseDir.resolve(args[0])
        if( !Files.isDirectory(moduleDir) )
            throw new AbortOperationException("Not a directory: ${moduleDir}")
        return moduleDir
    }

    private MetaYmlGenerator.ProcessMetadata parseModule(Path mainNfPath) {
        try {
            final metadata = MetaYmlGenerator.extract(mainNfPath)
            if( !metadata )
                throw new AbortOperationException("No process definition found in: ${mainNfPath}")
            return metadata
        }
        catch( AbortOperationException e ) {
            throw e
        }
        catch( Exception e ) {
            throw new AbortOperationException("Failed to parse ${mainNfPath}: ${e.message}")
        }
    }
}
