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
package nextflow.agent

import dev.langchain4j.model.chat.request.json.JsonArraySchema
import dev.langchain4j.model.chat.request.json.JsonBooleanSchema
import dev.langchain4j.model.chat.request.json.JsonIntegerSchema
import dev.langchain4j.model.chat.request.json.JsonNumberSchema
import dev.langchain4j.model.chat.request.json.JsonObjectSchema
import dev.langchain4j.model.chat.request.json.JsonSchema
import dev.langchain4j.model.chat.request.json.JsonSchemaElement
import dev.langchain4j.model.chat.request.json.JsonStringSchema
import groovy.transform.CompileStatic

/**
 * Converts a portable JSON-schema {@link Map} (the shape produced by core's
 * {@code nextflow.agent.RecordSchema.of}) into a langchain4j
 * {@link JsonSchema} suitable for use as a structured-output contract.
 *
 * The portable shape is:
 * <pre>
 * [ type:'object',
 *   properties:[ field:[ type:'string'|'integer'|'number'|'boolean'|'array'|'object',
 *                        items:..., properties:..., required:... ] ],
 *   required:[...],
 *   additionalProperties:false ]
 * </pre>
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class JsonSchemaMapper {

    /**
     * Build a langchain4j {@link JsonSchema} with the given name from a portable
     * object schema map.
     *
     * @param name   the schema name (used as the structured-output schema name)
     * @param schema the portable schema map (must describe an {@code object})
     * @return a non-null {@link JsonSchema} whose root is a {@link JsonObjectSchema}
     */
    static JsonSchema toJsonSchema(String name, Map schema) {
        final root = toObjectSchema(schema)
        return JsonSchema.builder()
            .name(name)
            .rootElement(root)
            .build()
    }

    private static JsonObjectSchema toObjectSchema(Map schema) {
        final builder = JsonObjectSchema.builder()
        final Map properties = (schema?.properties ?: [:]) as Map
        for( Map.Entry entry : properties.entrySet() ) {
            final key = entry.key as String
            final spec = entry.value as Map
            builder.addProperty(key, toElement(key, spec))
        }
        final required = schema?.required as List
        if( required != null )
            builder.required(required.collect { it as String })
        final additional = schema?.additionalProperties
        if( additional != null )
            builder.additionalProperties(additional as Boolean)
        return builder.build()
    }

    private static JsonSchemaElement toElement(String name, Map spec) {
        final type = spec?.type as String
        switch( type ) {
            case 'string':
                return new JsonStringSchema()
            case 'integer':
                return new JsonIntegerSchema()
            case 'number':
                return new JsonNumberSchema()
            case 'boolean':
                return new JsonBooleanSchema()
            case 'array':
                final items = spec.items as Map
                if( items == null )
                    throw new IllegalArgumentException("Array property `${name}` is missing an `items` schema")
                return JsonArraySchema.builder()
                    .items(toElement(name + '[]', items))
                    .build()
            case 'object':
                return toObjectSchema(spec)
            default:
                throw new IllegalArgumentException("Unsupported schema type `${type}` for property `${name}`")
        }
    }
}
