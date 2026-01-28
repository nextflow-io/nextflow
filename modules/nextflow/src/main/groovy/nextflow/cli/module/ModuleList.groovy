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
import groovy.json.JsonOutput
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.cli.CmdBase
import nextflow.exception.AbortOperationException
import nextflow.module.InstalledModule
import nextflow.module.ModuleIntegrity
import nextflow.module.ModuleStorage

import java.nio.file.Paths

/**
 * Module list subcommand
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Slf4j
@CompileStatic
@Parameters(commandDescription = "List all installed modules")
class ModuleList extends CmdBase {

    @Parameter(names = ["-json"], description = "Output in JSON format", arity=0)
    boolean jsonOutput = false

    @Override
    String getName() {
        return 'list'
    }

    @Override
    void run() {

        // Get config
        def baseDir = Paths.get('.').toAbsolutePath().normalize()


        // Create resolver and list modules
        def storage = new ModuleStorage(baseDir)

        try {
            def installed = storage.listInstalled()

            if (installed.isEmpty()) {
                println "No modules installed"
                return
            }

            if (jsonOutput) {
                printJsonList(installed)
            } else {
                printFormattedList(installed)
            }
        }
        catch (Exception e) {
            log.error("Failed to list modules", e)
            throw new AbortOperationException("List failed: ${e.message}", e)
        }
    }

    private void printFormattedList(List<InstalledModule> installed) {
        println ""
        println "Installed modules:"
        println ""
        println "Module".padRight(40) + "Version".padRight(15) + "Status"
        println ("-" * 70)

        installed.each { module ->
            def status = getStatusString(module.integrity)
            println "${module.reference.nameWithoutPrefix.padRight(40)}${(module.installedVersion ?: 'unknown').padRight(15)}${status}"
        }
        println ""
    }

    private void printJsonList(List<InstalledModule> installed) {
        def modules = installed.collect { module ->
            [
                name: module.reference.nameWithoutPrefix,
                version: module.installedVersion ?: 'unknown',
                integrity: module.integrity.toString(),
                directory: module.directory.toString()
            ]
        }

        // Simple JSON output (could use groovy.json.JsonOutput for better formatting)
        println JsonOutput.toJson(modules: modules)
    }

    private String getStatusString(ModuleIntegrity integrity) {
        switch (integrity) {
            case ModuleIntegrity.VALID:
                return 'OK'
            case ModuleIntegrity.MODIFIED:
                return 'MODIFIED'
            case ModuleIntegrity.MISSING_CHECKSUM:
                return 'NO CHECKSUM'
            case ModuleIntegrity.CORRUPTED:
                return 'CORRUPTED'
            default:
                return 'UNKNOWN'
        }
    }
}
