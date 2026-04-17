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

package nextflow.module

import java.nio.file.Files
import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Const

/**
 * Validates module structure and module spec (meta.yml).
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class ModuleValidator {

    /**
     * Run all validations on a module directory and return a list of error messages.
     * An empty list means the module is valid.
     *
     * @param moduleDir
     */
    static List<String> validate(Path moduleDir) {
        final errors = new ArrayList<String>()

        // Level 1: validate module structure
        errors.addAll(validateStructure(moduleDir))
        if( errors )
            return errors  // can't proceed without required files

        // Level 2: validate module spec (meta.yml)
        final manifestPath = moduleDir.resolve(ModuleStorage.MODULE_MANIFEST_FILE)
        final spec = ModuleSpecFactory.fromYaml(manifestPath)
        errors.addAll(spec.validate())
        if( errors )
            return errors

        // Level 3: validate module input/output spec against process definition
        final scriptPath = moduleDir.resolve("main.nf")
        final sourceSpec = ModuleSpecFactory.fromScript(scriptPath)
        errors.addAll(validateInputsOutputs(spec, sourceSpec))

        return errors
    }

    /**
     * Check that required files exist.
     *
     * @param moduleDir
     */
    static List<String> validateStructure(Path moduleDir) {
        final errors = new ArrayList<String>()

        if( !Files.exists(moduleDir) || !Files.isDirectory(moduleDir) ) {
            errors << "Module directory does not exist: ${moduleDir}".toString()
            return errors
        }

        if( !Files.exists(moduleDir.resolve(Const.DEFAULT_MAIN_FILE_NAME)) )
            errors << "Missing required file: ${Const.DEFAULT_MAIN_FILE_NAME}".toString()

        if( !Files.exists(moduleDir.resolve(ModuleStorage.MODULE_MANIFEST_FILE)) )
            errors << "Missing required file: ${ModuleStorage.MODULE_MANIFEST_FILE}".toString()

        if( !Files.exists(moduleDir.resolve(ModuleStorage.MODULE_README_FILE)) )
            errors << "Missing required file: ${ModuleStorage.MODULE_README_FILE}".toString()

        // Check bundle size (1MB uncompressed limit)
        try (final sizeStream = Files.walk(moduleDir)){
            long totalSize = sizeStream
                .filter { Files.isRegularFile(it) }
                .mapToLong { Files.size(it) }
                .sum()

            def maxSize = 1024 * 1024 // 1MB in bytes
            if (totalSize > maxSize) {
                def sizeMB = totalSize / (1024 * 1024)
                errors << "Module size exceeds 1MB limit (current: ${String.format('%.2f', sizeMB)}MB)".toString()
            }
        } catch (Exception e) {
            log.warn "Failed to check module size: ${e.message}"
        }

        return errors
    }

    static List<String> validateInputsOutputs(ModuleSpec spec, ModuleSpec sourceSpec) {
        final errors = new ArrayList<String>()

        final specInputs = spec.inputs?.size() ?: 0
        final sourceInputs = sourceSpec.inputs?.size() ?: 0
        if( specInputs != sourceInputs )
            errors << "Module spec has ${specInputs} inputs but module declares ${sourceInputs} inputs".toString()

        if( spec.outputs != null ) {
            final specOutputs = spec.outputs.size()
            final sourceOutputs = sourceSpec.outputs?.size() ?: 0
            if( specOutputs != sourceOutputs )
                errors << "Module spec has ${specOutputs} outputs but module declares ${sourceOutputs} outputs".toString()
        }

        if( spec.topics != null ) {
            final specTopics = spec.topics.size()
            final sourceTopics = sourceSpec.topics?.size() ?: 0
            if( specTopics != sourceTopics )
                errors << "Module spec has ${specTopics} topic emissions but module has ${sourceTopics} topic emissions".toString()
        }

        return errors
    }
}
