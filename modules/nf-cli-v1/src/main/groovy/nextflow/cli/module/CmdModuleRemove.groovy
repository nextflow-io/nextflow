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

import nextflow.util.TestOnly

import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths

/**
 * Module remove subcommand
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Slf4j
@CompileStatic
@Parameters(commandDescription = "Remove an installed module")
class CmdModuleRemove extends CmdBase {

    @Parameter(description = "<module>", required = true)
    List<String> args

    @Parameter(names = ["-keep-files"], description = "Remove only .module-info keeping the local files", arity = 0)
    boolean keepFiles = false

    @Parameter(names = ["-force"], description = "Force remove", arity = 0)
    boolean force = false

    @TestOnly
    protected Path root

    @Override
    String getName() {
        return 'remove'
    }

    @Override
    void run() {
        if( !args || args.size() != 1 ) {
            throw new AbortOperationException("Incorrect number of arguments")
        }
        if( keepFiles && force ) {
            throw new AbortOperationException("Cannot use both -keep-files and -force options")
        }

        def reference = ModuleReference.parse(args[0])

        // Get config
        def baseDir = root ?: Paths.get('.').toAbsolutePath().normalize()

        // Create resolver and spec file manager
        def storage = new ModuleStorage(baseDir)

        try {
            def filesRemoved = false

            // Remove local files unless -keep-files is set
            if( !keepFiles ) {
                filesRemoved = storage.removeModule(reference, force)
                if( filesRemoved ) {
                    println "Module ${reference} files removed successfully"
                } else {
                    println "Module ${reference} not found locally"
                }
            } else {
                println "Keeping module files for ${reference} (-keep-files flag)"
                final moduleInfo = storage.getModuleInfo(reference)
                if( Files.exists(moduleInfo) ) {
                    Files.delete(moduleInfo)
                }
            }

        }
        catch( AbortOperationException e ) {
            throw e
        }
        catch( Exception e ) {
            throw new AbortOperationException("Failed to remove module $reference: ${e.message}", e)
        }
    }
}
