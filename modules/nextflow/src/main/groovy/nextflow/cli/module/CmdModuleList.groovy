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

import com.beust.jcommander.IParameterValidator
import com.beust.jcommander.Parameter
import com.beust.jcommander.ParameterException
import com.beust.jcommander.Parameters
import groovy.json.JsonOutput
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.cli.CmdBase
import nextflow.exception.AbortOperationException
import nextflow.module.InstalledModule
import nextflow.module.ModuleIntegrity
import nextflow.module.ModuleStorage
import nextflow.util.TestOnly

import java.nio.file.Path
import java.nio.file.Paths

/**
 * Module list subcommand
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Slf4j
@CompileStatic
@Parameters(commandDescription = "List all installed modules")
class CmdModuleList extends CmdBase {

    @Parameter(
        names = ['-o', '-output'],
        description = 'Output mode for reporting search results: table, json',
        validateWith = OutputModeValidator
    )
    String output = 'table'

    static class OutputModeValidator implements IParameterValidator {

        private static final List<String> MODES = List.of('table', 'json')

        @Override
        void validate(String name, String value) {
            if( !MODES.contains(value) )
                throw new ParameterException("Output mode must be one of $MODES (found: $value)")
        }
    }

    @TestOnly
    protected Path root

    @Override
    String getName() {
        return 'list'
    }

    @Override
    void run() {

        // Get config
        def baseDir = root ?: Paths.get('.').toAbsolutePath().normalize()


        // Create resolver and list modules
        def storage = new ModuleStorage(baseDir)

        try {
            def installed = storage.listInstalled()

            if( installed.isEmpty() ) {
                println "No modules installed"
                return
            }

            if( !output || output == 'table' ) {
                printFormattedList(installed)
            } else if( output == 'json' ) {
                printJsonList(installed)
            } else {
                throw new AbortOperationException("Not implemented output mode $output)")
            }

        }
        catch( Exception e ) {
            log.error("Failed to list modules", e)
            throw new AbortOperationException("List failed: ${e.message}", e)
        }
    }

    private void printFormattedList(List<InstalledModule> installed) {
        println ""
        println "Installed modules:"
        println ""
        println "Module".padRight(40) + "Version".padRight(15) + "Status"
        println("-" * 70)

        installed.each { module ->
            def status = getStatusString(module.integrity)
            println "${module.reference.toString().padRight(40)}${(module.installedVersion ?: 'unknown').padRight(15)}${status}"
        }
        println ""
    }

    private void printJsonList(List<InstalledModule> installed) {
        def modules = installed.collect { module ->
            [
                name     : module.reference.toString(),
                version  : module.installedVersion ?: 'unknown',
                integrity: module.integrity.toString(),
                directory: module.directory.toString(),
                registry : module.registryUrl ?: 'unknown'
            ]
        }

        // Simple JSON output (could use groovy.json.JsonOutput for better formatting)
        println JsonOutput.prettyPrint(JsonOutput.toJson(modules: modules))
    }

    private String getStatusString(ModuleIntegrity integrity) {
        switch( integrity ) {
            case ModuleIntegrity.VALID:
                return 'OK'
            case ModuleIntegrity.MODIFIED:
                return 'MODIFIED'
            case ModuleIntegrity.NO_REMOTE_MODULE:
                return 'LOCAL'
            case ModuleIntegrity.CORRUPTED:
                return 'CORRUPTED'
            default:
                return 'UNKNOWN'
        }
    }
}
