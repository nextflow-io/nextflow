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
            def data = yaml.load(Files.newInputStream(metaYamlPath)) as Map<String, Object>

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
