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
import nextflow.module.ModuleStorage
import nextflow.pipeline.PipelineSpec

import java.nio.file.Paths

/**
 * Module remove subcommand
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Slf4j
@CompileStatic
@Parameters(commandDescription = "Remove an installed module")
class ModuleRemove extends CmdBase {

    @Parameter(description = "<module>", required = true)
    List<String> args

    @Parameter(names = ["-keep-config"], description = "Remove local files but keep the entry in nextflow_spec.json", arity = 0)
    boolean keepConfig = false

    @Parameter(names = ["-keep-files"], description = "Remove from config but keep local files", arity = 0)
    boolean keepFiles = false

    @Override
    String getName() {
        return 'remove'
    }

    @Override
    void run() {
        if (!args || args.size() != 1) {
            throw new AbortOperationException("Incorrect number of arguments")
        }

        // Validate flags
        if (keepConfig && keepFiles) {
            throw new AbortOperationException("Cannot use both -keep-config and -keep-files flags together")
        }

        def moduleRef = '@' + args[0]

        def reference = ModuleReference.parse(moduleRef)

        // Get config
        def baseDir = Paths.get('.').toAbsolutePath().normalize()

        //TODO: Decide final location of modules currently in nextflow_spec.json.
        def specFile = new PipelineSpec(baseDir)

        // Create resolver and spec file manager
        def storage = new ModuleStorage(baseDir)

        try {
            def filesRemoved = false
            def configRemoved = false

            // Remove local files unless -keep-files is set
            if (!keepFiles) {
                println "Removing module files for ${reference.nameWithoutPrefix}..."
                filesRemoved = storage.removeModule(reference)
                if (filesRemoved) {
                    println "Module files removed successfully"
                } else {
                    println "Module ${reference.nameWithoutPrefix} was not installed locally"
                }
            } else {
                println "Keeping module files for ${reference.nameWithoutPrefix} (due to -keep-files flag)"
            }

            // Remove config entry unless -keep-config is set
            if (!keepConfig) {
                println "Removing module entry from nextflow_spec.json..."
                configRemoved = specFile.removeModuleEntry(reference.fullName)
                if (configRemoved) {
                    println "Module entry removed from configuration"
                } else {
                    println "Module ${reference.nameWithoutPrefix} was not configured in nextflow_spec.json"
                }
            } else {
                println "Keeping module entry in nextflow_spec.json (due to -keep-config flag)"
            }

            // Summary
            if (filesRemoved || configRemoved) {
                println "\nModule ${reference.nameWithoutPrefix} removal completed"
            } else {
                println "\nModule ${reference.nameWithoutPrefix} was not found"
            }
        }
        catch (AbortOperationException e) {
            throw e
        }
        catch (Exception e) {
            log.error("Failed to remove module", e)
            throw new AbortOperationException("Removal failed: ${e.message}", e)
        }
    }
}
