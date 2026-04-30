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
import com.networknt.schema.ValidationMessage
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
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
        final schemaText = loadSchema(schemaLocation)
        final factory = JsonSchemaFactory.getInstance(SpecVersion.VersionFlag.V202012)
        final JsonSchema schema
        try {
            schema = factory.getSchema(schemaText)
        }
        catch( Exception e ) {
            throw new AbortOperationException("Invalid module schema at '${schemaLocation}': ${e.message}", e)
        }

        Object yamlData
        try( final stream = Files.newInputStream(metaYaml) ) {
            yamlData = new Yaml().load(stream)
        }
        catch( Exception e ) {
            throw new AbortOperationException("Failed to read module spec '${metaYaml}': ${e.message}", e)
        }

        final JsonNode node = JSON_MAPPER.valueToTree(yamlData)
        final Set<ValidationMessage> messages = schema.validate(node)
        return messages.collect { it.message }.toList()
    }

    static List<String> validate(Path metaYaml) {
        return validate(metaYaml, DEFAULT_SCHEMA_URL)
    }

    /**
     * Load the JSON schema text from a remote URL, file: URI, or local file path.
     * Hard-fails with AbortOperationException on any I/O error.
     */
    private static String loadSchema(String location) {
        try {
            if( location.startsWith('http://') || location.startsWith('https://') ) {
                final url = new URL(location)
                final conn = url.openConnection()
                conn.setConnectTimeout(10_000)
                conn.setReadTimeout(20_000)
                try( final stream = conn.getInputStream() ) {
                    return new String(stream.readAllBytes(), 'UTF-8')
                }
            }
            if( location.startsWith('file:') ) {
                return Files.readString(Paths.get(URI.create(location)))
            }
            return Files.readString(Paths.get(location))
        }
        catch( Exception e ) {
            throw new AbortOperationException(
                "Failed to load module schema from '${location}': ${e.message}. " +
                "Pass --schema <url-or-local-path> to override.", e)
        }
    }
}
