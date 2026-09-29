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

import dev.langchain4j.internal.JsonSchemaElementUtils
import dev.langchain4j.model.chat.request.json.JsonSchema
import spock.lang.Specification

/**
 * Characterization / regression guard for OpenAI-strict optional ({@code @Nullable}/{@code ?})
 * structured-output fields (milestone M4).
 *
 * The portable schema layer (core {@code RecordSchema.of}) omits {@code @Nullable} fields from
 * the {@code required} list. This test asserts what the OpenAI plugin actually sends on the wire:
 * langchain4j's strict serialization ({@code JsonSchemaElementUtils.toMap(root, true)} with
 * {@code strict=true}, exactly what {@code OpenAiUtils.toOpenAiResponseFormat} invokes) rewrites the object-level
 * {@code required} to include every property, BUT emits an OpenAI nullable-union type
 * ({@code ["<type>","null"]}) for fields omitted from {@code required}. That is the canonical
 * OpenAI-strict optional-field idiom: the field is listed as required but its type is nullable,
 * so the model may legitimately return it null. See {@code adr/specs/agent-design.md} §3.3.
 *
 * This lives in the plugin because the nullable-union is produced by langchain4j serialization,
 * which is on the plugin classpath only. It is a faithful, network-free proxy for the schema the
 * model receives. If a langchain4j version bump removes the {@code type()} nullable helper this
 * test turns red — the signal to apply the contingency source fix.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class NullableStrictSchemaTest extends Specification {

    /** Normalise a schema {@code type} which may be a plain {@code String} or a {@code String[]} union. */
    private static List asList(Object t) { t instanceof String[] ? (t as List) : (t instanceof List ? t : [t]) }

    private static Map toWire(String name, Map portable, boolean strict) {
        final JsonSchema schema = JsonSchemaMapper.toJsonSchema(name, portable)
        return JsonSchemaElementUtils.toMap(schema.rootElement(), strict)
    }

    def 'scalar optional field becomes a nullable union but stays in the strict required set'() {
        given: 'the portable map RecordSchema.of produces for `record Answer { answer: String; note: String? }`'
        def portable = [type                : 'object',
                        properties          : [answer: [type: 'string'], note: [type: 'string']],
                        required            : ['answer'],
                        additionalProperties: false]

        when: 'serialized with strict=true (as OpenAiUtils.toOpenAiResponseFormat does)'
        def wire = toWire('Answer', portable, true)

        then: 'the required field keeps a plain type'
        asList(wire.properties.answer.type) == ['string']

        and: 'the optional field carries the OpenAI nullable union'
        asList(wire.properties.note.type) == ['string', 'null']

        and: 'strict rewrote required to all keys - the field is "required" but nullable'
        (wire.required as Set) == ['answer', 'note'] as Set

        and:
        wire.additionalProperties == false
    }

    def 'optional field inside a nested record recurses to a nullable union'() {
        given: '`record Outer { title: String; inner: Inner }`, Inner { name: String; note: String? }'
        def portable = [type                : 'object',
                        properties          : [
                            title: [type: 'string'],
                            inner: [type                : 'object',
                                    properties          : [name: [type: 'string'], note: [type: 'string']],
                                    required            : ['name'],
                                    additionalProperties: false],
                        ],
                        required            : ['title', 'inner'],
                        additionalProperties: false]

        when:
        def wire = toWire('Outer', portable, true)

        then: 'the required nested object stays a plain object'
        asList(wire.properties.inner.type) == ['object']

        and: 'its optional field carries the nullable union at depth'
        asList(wire.properties.inner.properties.note.type) == ['string', 'null']

        and: 'the nested required set was rewritten to include the nullable field'
        (wire.properties.inner.required as Set) == ['name', 'note'] as Set
    }

    def 'a nullable nested object becomes an object/null union'() {
        given: '`record Outer { name: String; addr: Address }` with addr omitted from required (i.e. @Nullable)'
        def portable = [type                : 'object',
                        properties          : [
                            name: [type: 'string'],
                            addr: [type                : 'object',
                                   properties          : [city: [type: 'string']],
                                   required            : ['city'],
                                   additionalProperties: false],
                        ],
                        required            : ['name'],
                        additionalProperties: false]

        when:
        def wire = toWire('Outer', portable, true)

        then: 'the nullable object becomes an object/null union'
        asList(wire.properties.addr.type) == ['object', 'null']

        and: 'it is still listed in the strict required set'
        (wire.required as Set) == ['name', 'addr'] as Set
    }

    def 'the exit criterion holds at the M1 multi-output wrapper nesting level'() {
        given: 'a synthetic wrapper (independent of M1 code) whose record property carries a nullable field'
        def portable = [type                : 'object',
                        properties          : [
                            plan : [type                : 'object',
                                    properties          : [step: [type: 'string'], hint: [type: 'string']],
                                    required            : ['step'],
                                    additionalProperties: false],
                            count: [type: 'integer'],
                        ],
                        required            : ['plan', 'count'],
                        additionalProperties: false]

        when:
        def wire = toWire('Wrapper', portable, true)

        then: 'the deeply-nested nullable field still receives the nullable union - recursion at every wrapper level'
        asList(wire.properties.plan.properties.hint.type) == ['string', 'null']

        and:
        (wire.properties.plan.required as Set) == ['step', 'hint'] as Set
    }

    def 'non-strict serialization keeps the classic optional idiom (control)'() {
        given: 'the same single-record portable map'
        def portable = [type                : 'object',
                        properties          : [answer: [type: 'string'], note: [type: 'string']],
                        required            : ['answer'],
                        additionalProperties: false]

        when: 'serialized with strict=false'
        def wire = toWire('Answer', portable, false)

        then: 'the optional field keeps a plain scalar type'
        asList(wire.properties.note.type) == ['string']

        and: 'and is simply omitted from required (classic optional idiom)'
        !(wire.required as List).contains('note')
        (wire.required as List) == ['answer']
    }
}
