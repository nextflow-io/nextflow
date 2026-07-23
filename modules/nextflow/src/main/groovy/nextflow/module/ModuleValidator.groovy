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
     * @param schemaLocation URL or local path of the JSON schema used to validate meta.yml
     */
    static List<String> validate(Path moduleDir, String schemaLocation) {
        final errors = new ArrayList<String>()

        // Level 1: validate module structure
        errors.addAll(validateStructure(moduleDir))
        if( errors )
            return errors  // can't proceed without required files

        // Level 2a: validate module spec (meta.yml) against the JSON schema
        final manifestPath = moduleDir.resolve(ModuleStorage.MODULE_MANIFEST_FILE)
        errors.addAll(ModuleSchemaValidator.validate(manifestPath, schemaLocation))
        if( errors )
            return errors

        // Level 2b: validate Nextflow-specific rules not expressed by the schema
        final spec = ModuleSpecFactory.fromYaml(manifestPath)
        errors.addAll(spec.validate())
        if( errors )
            return errors

        // Level 2c: validate that declared module dependencies are vendored at their pinned version
        errors.addAll(validateRequiresModules(moduleDir, spec))
        if( errors )
            return errors

        // Level 3: validate the module script against its kind
        final scriptPath = moduleDir.resolve(Const.DEFAULT_MAIN_FILE_NAME)
        if( spec.isWorkflow() ) {
            errors.addAll(validateWorkflow(spec, scriptPath))
        }
        else {
            final sourceSpec = ModuleSpecFactory.fromScript(scriptPath)
            errors.addAll(validateInputsOutputs(spec, sourceSpec))
        }

        return errors
    }

    static List<String> validate(Path moduleDir) {
        final manifest = moduleDir.resolve(ModuleStorage.MODULE_MANIFEST_FILE)
        return validate(moduleDir, ModuleSchemaValidator.resolveSchemaLocation(manifest, null))
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

    /**
     * Validate that every {@code requires.modules} dependency declared in the meta.yml is
     * vendored under the module's own nested {@code modules/} directory at the pinned version.
     *
     * The publish bundle ships these vendored dependencies (nested per-module vendoring, ADR v2),
     * so a dependency that is declared but not vendored -- or vendored at a different version --
     * would produce a broken bundle. This check is local/offline (no registry access).
     *
     * @param moduleDir the module directory being validated
     * @param spec      the parsed module spec
     */
    static List<String> validateRequiresModules(Path moduleDir, ModuleSpec spec) {
        final errors = new ArrayList<String>()
        final deps = spec.requiresModules
        if( !deps )
            return errors

        final nested = new ModuleStorage(moduleDir)
        for( final dep : deps ) {
            final parsed = ModuleResolver.parseDependency(dep)
            final ref = parsed.reference
            final installed = nested.getInstalledModule(ref)
            if( installed == null ) {
                errors << ("Declared dependency '${dep}' is not vendored under modules/${ref.fullName} -- " +
                    "run `nextflow module install ${dep}` from the module directory").toString()
                continue
            }
            if( parsed.version && installed.installedVersion != parsed.version )
                errors << "Declared dependency '${dep}' is vendored at version ${installed.installedVersion} but ${parsed.version} is required".toString()
        }
        return errors
    }

    /**
     * Validate a workflow module's script and reconcile its {@code take:}/{@code emit:} interface
     * with the meta.yml. A workflow module must define exactly one workflow; when the meta.yml also
     * declares {@code input}/{@code output} (e.g. a typed workflow) the counts must match the
     * workflow's take/emit arity. When the meta.yml omits them the take:/emit: sections are the sole
     * source of truth and there is nothing to reconcile.
     *
     * @param spec       the parsed module spec
     * @param scriptPath the module's main.nf
     */
    static List<String> validateWorkflow(ModuleSpec spec, Path scriptPath) {
        final errors = new ArrayList<String>()

        final interfaces = ModuleSpecFactory.workflowInterfaces(scriptPath)
        if( interfaces.isEmpty() ) {
            errors << "Workflow module '${spec.name}' must define a workflow in ${Const.DEFAULT_MAIN_FILE_NAME}".toString()
            return errors
        }
        if( interfaces.size() > 1 ) {
            errors << "Workflow module '${spec.name}' must define exactly one workflow in ${Const.DEFAULT_MAIN_FILE_NAME} (found ${interfaces.size()})".toString()
            return errors
        }

        final iface = interfaces[0]
        if( spec.inputs != null && spec.inputs.size() != iface.takes )
            errors << "Module spec has ${spec.inputs.size()} inputs but workflow declares ${iface.takes} take(s)".toString()
        if( spec.outputs != null && spec.outputs.size() != iface.emits )
            errors << "Module spec has ${spec.outputs.size()} outputs but workflow declares ${iface.emits} emit(s)".toString()

        return errors
    }
}
