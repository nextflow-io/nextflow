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
import nextflow.exception.AbortOperationException
import nextflow.module.ModuleReference
import nextflow.module.ModuleSpec
import nextflow.module.ModuleStorage
import nextflow.util.TestOnly

import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths

/**
 * Implements the {@code nextflow module spec} subcommand.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Slf4j
@CompileStatic
@Parameters(commandDescription = "Generate a meta.yml manifest for a local module")
class CmdModuleSpec extends CmdBase {

    @Parameter(description = "Reference to local module (scope/name) or module directory path")
    List<String> args

    @Parameter(names = ["-force"], description = "Overwrite existing meta.yml", arity = 0)
    boolean force = false

    @Parameter(names = ["-dry-run"], description = "Print generated YAML to stdout without writing any file", arity = 0)
    boolean dryRun = false

    @Parameter(names = ["-scope"], description = "Module scope")
    String moduleScope

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

    private String moduleName

    @Override
    String getName() { 'spec' }

    @Override
    void run() {
        if( !args ){
            throw new AbortOperationException("Reference to local module (scope/name) or path is required")
        }
        final baseDir = root ?: Paths.get('.').toAbsolutePath().normalize()
        def moduleDir = resolveAsModuleReference(baseDir, args[0])
        if( moduleDir ) {
            moduleName = args[0]
            log.info("Module name inferred as $moduleName ${moduleScope ? '- scope parameter is ignored':''}")
        } else {
            moduleDir = resolveAsPath(baseDir, args[0])
            if( !moduleScope )
                throw new AbortOperationException("Module name can't be inferred from module directory path. Provide -scope parameter")
        }
        final mainNfPath = moduleDir.resolve('main.nf')

        if( !Files.exists(mainNfPath) )
            throw new AbortOperationException("Missing required file: main.nf in ${moduleDir}")

        final yamlContent = parseModule(mainNfPath).render()

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

    private Path resolveAsModuleReference(Path baseDir, String module) {
        try {
            final ref = ModuleReference.parse(module)
            final localStorage = new ModuleStorage(baseDir)
            if( localStorage.isInstalled(ref) ) {
                log.debug("Argument $module refers to a local module reference")
                return localStorage.getModuleDir(ref)
            }
            log.debug("Argument $module is not a local module reference")
            return null
        } catch( AbortOperationException e ) {
            log.debug("Argument $module is not a correct module reference")
            return null
        }
    }

    private Path resolveAsPath(Path baseDir, String module) {
        final path = Paths.get(module)
        final moduleDir = path.isAbsolute() ? path : baseDir.resolve(module)
        if( !Files.isDirectory(moduleDir) )
            throw new AbortOperationException("Not a directory: ${moduleDir}")
        return moduleDir
    }

    private ModuleSpec parseModule(Path mainNfPath) {
        try {
            final options = new ModuleSpec.ExtractOptions(
                name: moduleName,
                scope: moduleScope,
                version: moduleVersion,
                description: description,
                license: license,
                authors: authors
            )
            final metadata = ModuleSpec.extract(mainNfPath, options)
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
