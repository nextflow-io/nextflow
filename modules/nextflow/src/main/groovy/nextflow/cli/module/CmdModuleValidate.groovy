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
import java.nio.file.Paths

import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.cli.CmdBase
import nextflow.exception.AbortOperationException
import nextflow.module.ModuleReference
import nextflow.module.ModuleStorage
import nextflow.module.ModuleValidator
import nextflow.util.TestOnly

/**
 * Module validate subcommand - validates module structure and meta.yml consistency
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
@Parameters(commandDescription = "Validate a module structure and metadata")
class CmdModuleValidate extends CmdBase {

    @Parameter(description = "[namespace/name or path]", required = true)
    List<String> args

    @TestOnly
    protected Path root

    @Override
    String getName() {
        return 'validate'
    }

    @Override
    void run() {
        if( !args || args.size() != 1 )
            throw new AbortOperationException("Incorrect number of arguments - usage: nextflow module validate <namespace/name>")

        final moduleDir = determineModuleDir(args[0])
        final errors = ModuleValidator.validate(moduleDir)

        if( errors ) {
            throw new AbortOperationException(
                "Module validation failed:\n" + errors.collect { "  - ${it}" }.join('\n')
            )
        }

        println "Module validation passed."
    }

    /**
     * Resolve module path from argument (directory path or namespace/name reference)
     */
    protected Path determineModuleDir(String module) {
        final path = Paths.get(module)
        if( path.exists() )
            return path.toAbsolutePath().normalize()

        final ref = ModuleReference.parse(module)
        final storage = new ModuleStorage(root ?: Paths.get('.').toAbsolutePath().normalize())

        if( !storage.isInstalled(ref) )
            throw new AbortOperationException("No module directory found for $module")

        return storage.getModuleDir(ref)
    }
}
