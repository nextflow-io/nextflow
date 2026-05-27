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
import java.nio.file.Paths

import com.fasterxml.jackson.databind.JsonNode
import com.fasterxml.jackson.databind.ObjectMapper
import com.networknt.schema.JsonSchema
import com.networknt.schema.JsonSchemaFactory
import com.networknt.schema.SpecVersion
import com.networknt.schema.SpecVersionDetector
import com.networknt.schema.ValidationMessage
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.BuildInfo
import nextflow.exception.AbortOperationException
import org.yaml.snakeyaml.Yaml

/**
 * Validates a module spec (meta.yml) against the Nextflow module JSON schema.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class ModuleSchemaValidator {

    static final String DEFAULT_SCHEMA_URL =
        'https://raw.githubusercontent.com/nextflow-io/schemas/refs/heads/main/module/v1/schema.json'

    private static final ObjectMapper JSON_MAPPER = new ObjectMapper()

    /**
     * Validate a meta.yml file against the JSON schema located at the given
     * URL or local file path.
     *
     * @param metaYaml       Path to the meta.yml file to validate
     * @param schemaLocation URL (http/https), file: URI, or local file path of the schema
     * @return List of validation error messages, empty if the spec is valid
     */
    static List<String> validate(Path metaYaml, String schemaLocation) {
        final schemaNode = parseSchema(loadSchema(schemaLocation), schemaLocation)
        final specVersion = detectSpecVersion(schemaNode, schemaLocation)
        final schema = buildSchema(schemaNode, specVersion, schemaLocation)
        final metaNode = loadMeta(metaYaml)
        final Set<ValidationMessage> messages = schema.validate(metaNode)
        return messages.collect { it.message }.toList()
    }

    static List<String> validate(Path metaYaml) {
        return validate(metaYaml, DEFAULT_SCHEMA_URL)
    }

    private static JsonNode parseSchema(String schemaText, String schemaLocation) {
        try {
            return JSON_MAPPER.readTree(schemaText)
        }
        catch( Exception e ) {
            throw new AbortOperationException("Invalid module schema at '${schemaLocation}': ${e.message}", e)
        }
    }

    private static SpecVersion.VersionFlag detectSpecVersion(JsonNode schemaNode, String schemaLocation) {
        try {
            return SpecVersionDetector.detect(schemaNode)
        }
        catch( Exception e ) {
            throw new AbortOperationException(
                "Cannot determine JSON Schema draft for '${schemaLocation}': ${e.message}. " +
                "The schema must declare a supported \$schema (e.g. https://json-schema.org/draft/2020-12/schema).", e)
        }
    }

    private static JsonSchema buildSchema(JsonNode schemaNode, SpecVersion.VersionFlag specVersion, String schemaLocation) {
        try {
            return JsonSchemaFactory.getInstance(specVersion).getSchema(schemaNode)
        }
        catch( Exception e ) {
            throw new AbortOperationException("Invalid module schema at '${schemaLocation}': ${e.message}", e)
        }
    }

    private static JsonNode loadMeta(Path metaYaml) {
        try( final stream = Files.newInputStream(metaYaml) ) {
            return JSON_MAPPER.valueToTree(new Yaml().load(stream))
        }
        catch( Exception e ) {
            throw new AbortOperationException("Failed to read module spec '${metaYaml}': ${e.message}", e)
        }
    }

    /**
     * Load the JSON schema text from a remote URL, file: URI, or local file path.
     * Hard-fails with AbortOperationException on any I/O error.
     */
    private static String loadSchema(String location) {
        final lower = location.toLowerCase()
        try {
            if( lower.startsWith('http://') || lower.startsWith('https://') ) {
                final url = new URL(location)
                final conn = url.openConnection()
                conn.setConnectTimeout(10_000)
                conn.setReadTimeout(20_000)
                conn.setRequestProperty('User-Agent', "Nextflow/${BuildInfo.version}")
                try( final stream = conn.getInputStream() ) {
                    return new String(stream.readAllBytes(), 'UTF-8')
                }
            }
            if( lower.startsWith('file:') ) {
                return Files.readString(Paths.get(URI.create(location)))
            }
            return Files.readString(Paths.get(location))
        }
        catch( Exception e ) {
            throw new AbortOperationException(
                "Failed to load module schema from '${location}': ${e.message}. " +
                "Pass -schema <url-or-local-path> to override.", e)
        }
    }
}
