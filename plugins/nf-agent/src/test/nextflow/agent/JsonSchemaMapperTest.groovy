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
import dev.langchain4j.model.chat.request.json.JsonStringSchema
import spock.lang.Specification

class JsonSchemaMapperTest extends Specification {

    def 'should map a flat object schema to a JsonSchema'() {
        given:
        def schema = [
            type: 'object',
            properties: [
                answer    : [type: 'string'],
                confidence: [type: 'number'],
            ],
            required: ['answer'],
            additionalProperties: false,
        ]

        when:
        JsonSchema result = JsonSchemaMapper.toJsonSchema('Answer', schema)

        then:
        result != null
        result.name() == 'Answer'
        and:
        result.rootElement() instanceof JsonObjectSchema
        def root = result.rootElement() as JsonObjectSchema
        root.properties().keySet() == ['answer', 'confidence'] as Set
        root.properties().get('answer') instanceof JsonStringSchema
        root.properties().get('confidence') instanceof JsonNumberSchema
        root.required() == ['answer']
        root.additionalProperties() == Boolean.FALSE
    }

    def 'should map all scalar property types'() {
        given:
        def schema = [
            type: 'object',
            properties: [
                s: [type: 'string'],
                i: [type: 'integer'],
                n: [type: 'number'],
                b: [type: 'boolean'],
            ],
            required: ['s', 'i', 'n', 'b'],
        ]

        when:
        def root = JsonSchemaMapper.toJsonSchema('All', schema).rootElement() as JsonObjectSchema

        then:
        root.properties().get('s') instanceof JsonStringSchema
        root.properties().get('i') instanceof JsonIntegerSchema
        root.properties().get('n') instanceof JsonNumberSchema
        root.properties().get('b') instanceof JsonBooleanSchema
    }

    def 'should map an array property recursing on items'() {
        given:
        def schema = [
            type: 'object',
            properties: [
                tags: [type: 'array', items: [type: 'string']],
            ],
        ]

        when:
        def root = JsonSchemaMapper.toJsonSchema('Tagged', schema).rootElement() as JsonObjectSchema

        then:
        root.properties().get('tags') instanceof JsonArraySchema
        (root.properties().get('tags') as JsonArraySchema).items() instanceof JsonStringSchema
    }

    def 'should map a nested object property recursively'() {
        given:
        def schema = [
            type: 'object',
            properties: [
                inner: [
                    type: 'object',
                    properties: [name: [type: 'string']],
                    required: ['name'],
                ],
            ],
        ]

        when:
        def root = JsonSchemaMapper.toJsonSchema('Outer', schema).rootElement() as JsonObjectSchema

        then:
        root.properties().get('inner') instanceof JsonObjectSchema
        def inner = root.properties().get('inner') as JsonObjectSchema
        inner.properties().get('name') instanceof JsonStringSchema
        inner.required() == ['name']
    }

    def 'should fail on an unsupported property type'() {
        when:
        JsonSchemaMapper.toJsonSchema('Bad', [type: 'object', properties: [x: [type: 'whatever']]])

        then:
        thrown(IllegalArgumentException)
    }
}
