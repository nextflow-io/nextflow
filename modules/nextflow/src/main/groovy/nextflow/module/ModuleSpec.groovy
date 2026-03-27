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

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.exception.AbortOperationException
import org.yaml.snakeyaml.Yaml

import java.nio.file.Files
import java.nio.file.Path

/**
 * Represents a module spec (meta.yml) with validation
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Slf4j
@CompileStatic
class ModuleSpec {

    String name
    String version
    String description
    List<String> authors
    String license
    List<String> keywords
    Map<String, String> requires

    /**
     * Load a module spec from a file
     *
     * @param path Path to module spec file
     * @return ModuleSpec instance
     */
    static ModuleSpec load(Path path) {
        if( !Files.exists(path) ) {
            throw new AbortOperationException("Module spec not found: ${path}")
        }

        try {
            def yaml = new Yaml()
            Map<String, Object> data
            try( def stream = Files.newInputStream(path) ) {
                data = yaml.load(stream) as Map<String, Object>
            }

            def spec = new ModuleSpec()
            spec.name = data.name as String
            spec.version = data.version as String
            spec.description = data.description as String
            spec.authors = data.authors as List<String> ?: []
            spec.license = data.license as String
            spec.keywords = data.keywords as List<String> ?: []
            spec.requires = data.requires as Map<String, String> ?: [:]

            return spec
        }
        catch( Exception e ) {
            throw new AbortOperationException("Failed to parse module spec: ${path}", e)
        }
    }

    /**
     * Load a map of input name -> declared type from a module spec.
     *
     * @param path Path to module spec file
     * @return Map of input name -> type, or null if the file has no input section
     */
    static Map<String, Class> loadInputTypes(Path path) {
        if( !Files.exists(path) )
            return Collections.emptyMap()

        Map<String, Object> data
        try( final stream = Files.newInputStream(path) ) {
            data = new Yaml().load(stream) as Map<String, Object>
        }
        catch( Exception e ) {
            log.warn "Failed to parse module spec at ${path}: ${e.message}"
            return Collections.emptyMap()
        }

        final inputSection = data?.get('input') as List
        if( !inputSection )
            return Collections.emptyMap()

        final Map<String, Class> inputTypes = new LinkedHashMap<>()
        for( final item : inputSection )
            extractInputTypes(item, inputTypes)

        return inputTypes
    }

    /**
     * Recursively extract name -> type entries from a structuredParameter item.
     *
     * Handles:
     * - Tuples: recursively traverse tuple components
     * - New module spec format: {name: "...", type: "..."}
     * - Old nf-core format: {<name>: {type: "...", ...}}
     *
     * @param item
     * @param result
     */
    private static void extractInputTypes(Object item, Map<String, Class> result) {
        if( item instanceof List ) {
            // Traverse tuple components as individual params
            for( final element : (List) item ) {
                extractInputTypes(element, result)
            }
        }
        else if( item instanceof Map ) {
            final m = (Map) item
            if( m.values().size() == 1 && m.values().first() instanceof Map ) {
                // Old nf-core format: {<name>: {type: "...", description: "..."}}
                final name = m.keySet().first()
                final type = ((Map) m[name]).get('type')
                if( type != null )
                    result.put(name.toString(), inputType(type.toString()))
            }
            else {
                // New paramSpec format: {name: "...", type: "..."}
                final name = m.get('name')
                final type = m.get('type')
                if( name != null && type != null )
                    result.put(name.toString(), inputType(type.toString()))
            }
        }
    }

    /**
     * Return the corresponding class for a given type name.
     *
     * @param type
     */
    private static Class inputType(String type) {
        return switch( type ) {
            case 'boolean' -> Boolean
            case 'file' -> Path
            case 'directory' -> Path
            case 'float' -> Number
            case 'integer' -> Integer
            case 'map' -> Map
            case 'string' -> String
            default -> null
        }
    }

    /**
     * Validate the module spec for required fields
     *
     * @return List of validation errors (empty if valid)
     */
    List<String> validate() {
        List<String> errors = []

        if( !name ) {
            errors << "Missing required field: name"
        }
        if( !version ) {
            errors << "Missing required field: version"
        }
        if( !description ) {
            errors << "Missing required field: description"
        }
        if( !license ) {
            errors << "Missing required field: license"
        }

        // Validate version format (semantic versioning)
        if( version && !version.matches(/^\d+\.\d+\.\d+(-[\w.-]+)?$/) ) {
            errors << "Invalid version format: ${version} (expected semantic versioning, e.g., 1.0.0)".toString()
        }

        // Validate name format (scope/name or scope/path/to/name for nested modules)
        if( name && !name.matches(/^[a-zA-Z0-9._-]+\/[a-zA-Z0-9_-]+(?:\/[a-zA-Z0-9_-]+)*$/) ) {
            errors << "Invalid module name format: ${name} (expected scope/name or scope/path/to/name, e.g., nf-core/fastqc or nf-core/gfatools/gfa2fa)".toString()
        }

        return errors
    }

    /**
     * Check if the spec is valid
     *
     * @return true if valid, false otherwise
     */
    boolean isValid() {
        return validate().isEmpty()
    }
}
