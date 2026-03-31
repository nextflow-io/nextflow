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
import java.util.regex.Pattern

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Const

/**
 * Validates module structure, metadata schema, and cross-validates
 * meta.yml input declarations against the process definition.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class ModuleValidator {

    /**
     * Run all validations on a module directory and return a list of error messages.
     * An empty list means the module is valid.
     */
    static List<String> validate(Path moduleDir) {
        final errors = new ArrayList<String>()

        // Level 1: structural validation
        errors.addAll(validateStructure(moduleDir))
        if( errors )
            return errors  // can't proceed without required files

        // Level 2: meta.yml schema validation
        final manifestPath = moduleDir.resolve(ModuleStorage.MODULE_MANIFEST_FILE)
        final spec = ModuleSpec.load(manifestPath)
        errors.addAll(spec.validate())
        if( errors )
            return errors  // can't cross-validate with a malformed spec

        // Level 3: cross-validate meta.yml inputs with process definition
        final mainNf = moduleDir.resolve(Const.DEFAULT_MAIN_FILE_NAME)
        errors.addAll(crossValidateInputs(mainNf, manifestPath))

        return errors
    }

    /**
     * Check that required files exist.
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

    /**
     * Cross-validate meta.yml input declarations against the process definition.
     */
    static List<String> crossValidateInputs(Path mainNf, Path manifestPath) {
        final errors = new ArrayList<String>()

        // load declared inputs from meta.yml
        final metaInputs = ModuleSpec.loadInputTypes(manifestPath)
        if( !metaInputs )
            return errors  // no inputs declared — nothing to cross-validate

        // extract process input names from script text
        final processInputNames = extractProcessInputNames(mainNf.text)
        if( !processInputNames )
            return errors  // no process found or no inputs — skip

        // check for inputs declared in meta.yml but not in the process
        for( final name : metaInputs.keySet() ) {
            if( !processInputNames.contains(name) )
                errors << "Input '${name}' declared in meta.yml but not found in process definition".toString()
        }

        // check for inputs in the process but not declared in meta.yml
        for( final name : processInputNames ) {
            if( !metaInputs.containsKey(name) )
                errors << "Input '${name}' found in process definition but not declared in meta.yml".toString()
        }

        return errors
    }

    /**
     * Extract process input parameter names from script text.
     *
     * Parses input declarations such as:
     *   val name, path name, val(name), path(name),
     *   tuple val(name1), path(name2), env name, each name
     */
    static Set<String> extractProcessInputNames(String scriptText) {
        final names = new LinkedHashSet<String>()

        // find the input: block inside a process definition
        final inputBlockPattern = Pattern.compile(
            /(?s)process\s+\w+\s*\{.*?input:\s*\n(.*?)(?=\n\s*(?:output|script|exec|shell|stub|when):)/
        )
        final matcher = inputBlockPattern.matcher(scriptText)
        if( !matcher.find() )
            return names

        final inputBlock = matcher.group(1)

        // extract parameter names from input declarations
        final paramPattern = Pattern.compile(
            /(?:val|path|env|each)\s*\(?\s*(\w+)\s*\)?/
        )
        final paramMatcher = paramPattern.matcher(inputBlock)
        while( paramMatcher.find() ) {
            names.add(paramMatcher.group(1))
        }

        return names
    }
}
