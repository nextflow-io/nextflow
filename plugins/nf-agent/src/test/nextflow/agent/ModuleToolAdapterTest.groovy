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

import dev.langchain4j.model.chat.request.json.JsonObjectSchema
import dev.langchain4j.model.chat.request.json.JsonStringSchema
import spock.lang.Specification

class ModuleToolAdapterTest extends Specification {

    def 'should build a ToolSpecification with a JsonObjectSchema param schema'() {
        given:
        def descriptor = new ToolDescriptor(
            'greet',
            'greet someone by name',
            [
                type: 'object',
                properties: [name: [type: 'string']],
                required: ['name'],
                additionalProperties: false,
            ],
            [type: 'object', properties: [answer: [type: 'string']]])

        when:
        def spec = ModuleToolAdapter.toToolSpecification(descriptor)

        then:
        spec.name() == 'greet'
        spec.description() == 'greet someone by name'

        and: 'the params are a JsonObjectSchema exposing the declared property'
        spec.parameters() instanceof JsonObjectSchema
        def params = spec.parameters() as JsonObjectSchema
        params.properties().containsKey('name')
        params.properties().get('name') instanceof JsonStringSchema
        params.required() == ['name']
        params.additionalProperties() == false
    }

    def 'should tolerate a descriptor with no input properties'() {
        given:
        def descriptor = new ToolDescriptor('noop', 'does nothing', [type: 'object'], null)

        when:
        def spec = ModuleToolAdapter.toToolSpecification(descriptor)

        then:
        spec.name() == 'noop'
        spec.parameters() instanceof JsonObjectSchema
    }

    def 'should fail when the descriptor is null'() {
        when:
        ModuleToolAdapter.toToolSpecification(null)
        then:
        thrown(IllegalArgumentException)
    }

    def 'should fail when the descriptor name is missing'() {
        when:
        ModuleToolAdapter.toToolSpecification(new ToolDescriptor(null, 'd', [type: 'object'], null))
        then:
        thrown(IllegalArgumentException)
    }
}
