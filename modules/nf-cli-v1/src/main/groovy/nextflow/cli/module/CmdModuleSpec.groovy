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

import java.nio.file.Files
import java.nio.file.Path

import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.cli.CmdBase
import nextflow.exception.AbortOperationException
import nextflow.module.ModuleReference
import nextflow.module.ModuleSpec
import nextflow.module.ModuleSpecFactory
import nextflow.module.ModuleStorage
import nextflow.util.TestOnly

/**
 * Implements the {@code nextflow module spec} subcommand.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Slf4j
@CompileStatic
@Parameters(commandDescription = "Generate a meta.yml spec for a local module")
class CmdModuleSpec extends CmdBase {

    @Parameter(description = "Reference to local module (namespace/name) or module directory path")
    List<String> args

    @Parameter(names = ["-dry-run"], description = "Print generated YAML to stdout without writing any file", arity = 0)
    boolean dryRun = false

    @Parameter(names = ["-namespace"], description = "Module namespace")
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
            throw new AbortOperationException("Reference to local module (namespace/name) or path is required")
        }
        final baseDir = root ?: Path.of('.').toAbsolutePath().normalize()
        def moduleDir = resolveAsModuleReference(baseDir, args[0])
        if( moduleDir ) {
            moduleName = args[0]
            if( !dryRun )
                log.info "Module name inferred as ${moduleName} ${moduleScope ? '(-namespace option is ignored)' : ''}"
        }
        else {
            moduleDir = resolveAsPath(baseDir, args[0])
            if( !moduleScope )
                throw new AbortOperationException("The -namespace option is required when referencing a module by path")
        }

        final mainNfPath = moduleDir.resolve('main.nf')
        if( !Files.exists(mainNfPath) )
            throw new AbortOperationException("Missing module script (main.nf) in ${moduleDir}")

        final metaYmlPath = moduleDir.resolve('meta.yml')
        if( Files.exists(metaYmlPath) && !dryRun )
            log.info "Using existing module spec: ${metaYmlPath}"

        final yamlContent = generateSpec(mainNfPath, metaYmlPath).toYaml()

        if( dryRun ) {
            println yamlContent
            return
        }

        Files.writeString(metaYmlPath, yamlContent)
        println "Saved module spec to: ${metaYmlPath}"
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
        final path = Path.of(module)
        final moduleDir = path.isAbsolute() ? path : baseDir.resolve(module)
        if( !Files.isDirectory(moduleDir) )
            throw new AbortOperationException("Not a directory: ${moduleDir}")
        return moduleDir
    }

    private ModuleSpec generateSpec(Path mainNfPath, Path metaYmlPath) {
        final oldSpec = Files.exists(metaYmlPath)
            ? ModuleSpecFactory.fromYaml(metaYmlPath)
            : new ModuleSpec()
        return ModuleSpecFactory.fromScript(
            mainNfPath,
            oldSpec,
            namespace: moduleScope,
            name: moduleName,
            version: moduleVersion,
            description: description,
            license: license,
            authors: authors
        )
    }
}
