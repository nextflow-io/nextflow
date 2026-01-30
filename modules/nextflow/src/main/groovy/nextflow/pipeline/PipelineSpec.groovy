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

package nextflow.pipeline

import groovy.json.JsonOutput
import groovy.json.JsonSlurper
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j

import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.StandardOpenOption

/**
 * Manages the nextflow_spec.json file for module version declarations
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Slf4j
@CompileStatic
class PipelineSpec {

    private static final String SPEC_FILE_NAME = 'nextflow_spec.json'

    private final Path baseDir
    private final Path specFile

    PipelineSpec(Path baseDir) {
        this.baseDir = baseDir
        this.specFile = baseDir.resolve(SPEC_FILE_NAME)
    }

    /**
     * Add or update a module entry in the spec file
     *
     * @param moduleName The module name (e.g., "@nf-core/fastqc" or "nf-core/fastqc")
     * @param version The module version
     */
    void addModuleEntry(String moduleName, String version) {
        // Normalize module name (ensure it starts with @)
        def normalizedName = moduleName.startsWith('@') ? moduleName : '@' + moduleName

        def spec = readSpecFile()

        // Ensure modules map exists
        if (!spec.modules) {
            spec.modules = [:]
        }

        // Check if already configured with same version
        if (spec.modules[normalizedName] == version) {
            log.info "Module ${normalizedName} already configured with version ${version} in ${SPEC_FILE_NAME}"
            return
        }

        // Add or update entry
        spec.modules[normalizedName] = version
        writeSpecFile(spec)
        log.info "Added ${normalizedName}@${version} to ${SPEC_FILE_NAME}"
    }

    /**
     * Remove a module entry from the spec file
     *
     * @param moduleName The module name (e.g., "@nf-core/fastqc" or "nf-core/fastqc")
     * @return true if entry was removed, false if it didn't exist
     */
    boolean removeModuleEntry(String moduleName) {
        // Normalize module name (ensure it starts with @)
        def normalizedName = moduleName.startsWith('@') ? moduleName : '@' + moduleName

        def spec = readSpecFile()

        if (!spec.modules) {
             return false
        }
        final modules = spec.modules as Map<String,String>
        modules.remove(normalizedName)
        writeSpecFile(spec)
        log.info "Removed ${normalizedName} from ${SPEC_FILE_NAME}"
        return true
    }
    /**
     * @return Modules Map stored in the spec file
     */
    Map<String,String> getModules() {
        def spec = readSpecFile()

        if (!spec.modules) {
             return [:]
        }
        return spec.modules as Map<String,String>
    }

    /**
     * Check if the spec file exists
     *
     * @return true if the file exists
     */
    boolean exists() {
        return Files.exists(specFile)
    }

    private Map readSpecFile() {
        if (!Files.exists(specFile)) {
            return [:]
        }

        try {
            def content = Files.readString(specFile)
            if (content.trim().isEmpty()) {
                return [:]
            }
            return new JsonSlurper().parseText(content) as Map
        } catch (Exception e) {
            throw new RuntimeException("Failed to read spec file ${specFile}: ${e.message}", e)
        }
    }

    private void writeSpecFile(Map spec) {
        try {
            // Create directory if it doesn't exist
            if (!Files.exists(specFile.parent)) {
                Files.createDirectories(specFile.parent)
            }

            def jsonContent = JsonOutput.prettyPrint(JsonOutput.toJson(spec))
            Files.writeString(specFile, jsonContent, StandardOpenOption.CREATE, StandardOpenOption.TRUNCATE_EXISTING)
        } catch (Exception e) {
            throw new RuntimeException("Failed to write spec file ${specFile}: ${e.message}", e)
        }
    }
}