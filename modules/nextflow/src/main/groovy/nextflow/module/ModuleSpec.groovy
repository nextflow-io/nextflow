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
 * Represents a module manifest (meta.yaml) with validation
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
     * Load a module manifest from a meta.yaml file
     *
     * @param metaYamlPath Path to meta.yaml
     * @return ModuleSpec instance
     */
    static ModuleSpec load(Path metaYamlPath) {
        if (!Files.exists(metaYamlPath)) {
            throw new AbortOperationException("Module manifest not found: ${metaYamlPath}")
        }

        try {
            def yaml = new Yaml()
            def data = yaml.load(Files.newInputStream(metaYamlPath)) as Map<String, Object>

            def manifest = new ModuleSpec()
            manifest.name = data.name as String
            manifest.version = data.version as String
            manifest.description = data.description as String
            manifest.authors = data.authors as List<String> ?: []
            manifest.license = data.license as String
            manifest.keywords = data.keywords as List<String> ?: []
            manifest.requires = data.requires as Map<String, String> ?: [:]

            return manifest
        }
        catch (Exception e) {
            throw new AbortOperationException("Failed to parse module manifest: ${metaYamlPath}", e)
        }
    }

    /**
     * Validate the module manifest for required fields
     *
     * @return List of validation errors (empty if valid)
     */
    List<String> validate() {
        List<String> errors = []

        if (!name) {
            errors << "Missing required field: name"
        }
        if (!version) {
            errors << "Missing required field: version"
        }
        if (!description) {
            errors << "Missing required field: description"
        }
        if (!license) {
            errors << "Missing required field: license"
        }

        // Validate version format (semantic versioning)
        if (version && !version.matches(/^\d+\.\d+\.\d+(-[\w.-]+)?$/)) {
            errors << "Invalid version format: ${version} (expected semantic versioning, e.g., 1.0.0)".toString()
        }

        // Validate name format (scope/name)
        if (name && !name.matches(/^[a-zA-Z0-9_-]+\/[a-zA-Z0-9_-]+$/)) {
            errors << "Invalid module name format: ${name} (expected scope/name, e.g., nf-core/fastqc)".toString()
        }

        return errors
    }

    /**
     * Check if the manifest is valid
     *
     * @return true if valid, false otherwise
     */
    boolean isValid() {
        return validate().isEmpty()
    }
}
