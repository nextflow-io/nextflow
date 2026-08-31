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
import dev.langchain4j.model.chat.request.json.JsonEnumSchema
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

    def 'should propagate per-field descriptions, nested object properties and enum'() {
        given:
        def schema = [
            type: 'object',
            properties: [
                meta : [
                    type                : 'object',
                    description         : 'Groovy Map ... e.g. [id:..]',
                    properties          : [id: [type: 'string', description: 'sample identifier']],
                    additionalProperties: true,
                ],
                reads: [type: 'string', description: 'input reads (file path)'],
            ],
            required: ['meta', 'reads'],
            additionalProperties: false,
        ]

        when:
        def root = JsonSchemaMapper.toObjectSchema(schema)

        then: 'scalar property carries its description'
        def reads = root.properties().get('reads') as JsonStringSchema
        reads.description() == 'input reads (file path)'

        and: 'object property recurses and carries its description'
        def meta = root.properties().get('meta') as JsonObjectSchema
        meta instanceof JsonObjectSchema
        meta.description() == 'Groovy Map ... e.g. [id:..]'
        meta.additionalProperties() == Boolean.TRUE
        def id = meta.properties().get('id') as JsonStringSchema
        id instanceof JsonStringSchema
        id.description() == 'sample identifier'

        and: 'root structure preserved'
        root.required() == ['meta', 'reads']
        root.additionalProperties() == Boolean.FALSE
    }

    def 'should map a string with enum to an enum schema'() {
        given:
        def schema = [
            type: 'object',
            properties: [
                mode: [type: 'string', enum: ['fast', 'slow'], description: 'run mode'],
            ],
        ]

        when:
        def root = JsonSchemaMapper.toObjectSchema(schema)

        then:
        def mode = root.properties().get('mode')
        mode instanceof JsonEnumSchema
        (mode as JsonEnumSchema).enumValues() == ['fast', 'slow']
        (mode as JsonEnumSchema).description() == 'run mode'
    }

    def 'should propagate descriptions on all scalar types and arrays'() {
        given:
        def schema = [
            type: 'object',
            properties: [
                s   : [type: 'string', description: 'a string'],
                i   : [type: 'integer', description: 'an integer'],
                n   : [type: 'number', description: 'a number'],
                b   : [type: 'boolean', description: 'a boolean'],
                tags: [type: 'array', description: 'a list', items: [type: 'string']],
            ],
        ]

        when:
        def root = JsonSchemaMapper.toObjectSchema(schema)

        then:
        (root.properties().get('s') as JsonStringSchema).description() == 'a string'
        (root.properties().get('i') as JsonIntegerSchema).description() == 'an integer'
        (root.properties().get('n') as JsonNumberSchema).description() == 'a number'
        (root.properties().get('b') as JsonBooleanSchema).description() == 'a boolean'
        (root.properties().get('tags') as JsonArraySchema).description() == 'a list'
        (root.properties().get('tags') as JsonArraySchema).items() instanceof JsonStringSchema
    }
}
