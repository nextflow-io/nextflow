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
     * Load a module spec from a meta.yml file
     *
     * @param metaYamlPath Path to meta.yml
     * @return ModuleSpec instance
     */
    static ModuleSpec load(Path metaYamlPath) {
        if( !Files.exists(metaYamlPath) ) {
            throw new AbortOperationException("Module spec not found: ${metaYamlPath}")
        }

        try {
            def yaml = new Yaml()
            Map<String, Object> data
            try( def stream = Files.newInputStream(metaYamlPath) ) {
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
            throw new AbortOperationException("Failed to parse module spec: ${metaYamlPath}", e)
        }
    }

    /**
     * Load a flat map of input parameter name -> declared type from a meta.yml file.
     * Both the paramSpec format ({name, type, description}) and the old nf-core format
     * ({paramName: {type, description}}) are supported. Tuple inputs are flattened.
     *
     * @param metaYamlPath Path to meta.yml
     * @return Map of parameter name -> type string, or null if the file has no input section
     */
    static Map<String, String> loadInputTypes(Path metaYamlPath) {
        if( !Files.exists(metaYamlPath) ) return null

        Map<String, Object> data
        try {
            try( final stream = Files.newInputStream(metaYamlPath) ) {
                data = new Yaml().load(stream) as Map<String, Object>
            }
        }
        catch( Exception e ) {
            log.warn "Failed to parse meta.yml at ${metaYamlPath}: ${e.message}"
            return null
        }

        final inputSection = data?.get('input') as List
        if( !inputSection ) return null

        final Map<String, String> typeMap = new LinkedHashMap<>()
        for( final item : inputSection ) {
            extractParamTypes(item, typeMap)
        }

        return typeMap.isEmpty() ? null : typeMap
    }

    /**
     * Recursively extracts parameter name -> type entries from a structuredParameter item.
     * Handles:
     * - paramSpec format: {name: "...", type: "..."}
     * - Old nf-core format: {paramName: {type: "...", ...}}
     * - Arrays (tuples): recursively flatten all items
     */
    private static void extractParamTypes(Object item, Map<String, String> result) {
        if( item instanceof List ) {
            for( final element : (List) item ) {
                extractParamTypes(element, result)
            }
        }
        else if( item instanceof Map ) {
            final m = (Map) item
            final name = m.get('name')
            final type = m.get('type')
            if( name != null && type != null ) {
                // New paramSpec format: {name: "...", type: "..."}
                result.put(name.toString(), type.toString())
            }
            else {
                // Old nf-core format: {paramName: {type: "...", description: "..."}}
                // Also handles tuple-as-map-value: {tuple: [{paramName: {type: ...}}, ...]}
                for( final entry : ((Map<Object, Object>) m).entrySet() ) {
                    if( entry.value instanceof Map ) {
                        final innerType = ((Map) entry.value).get('type')
                        if( innerType != null ) {
                            result.put(entry.key.toString(), innerType.toString())
                        }
                    }
                    else if( entry.value instanceof List ) {
                        // Tuple described as a map value, e.g. {tuple: [{meta: {type: map}}, ...]}
                        for( final element : (List) entry.value ) {
                            extractParamTypes(element, result)
                        }
                    }
                }
            }
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
