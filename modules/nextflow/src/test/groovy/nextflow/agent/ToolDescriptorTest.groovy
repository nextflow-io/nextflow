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

import spock.lang.Specification

/**
 * Verifies the langchain4j-free agent tool SPI: the {@link ToolDescriptor} DTO
 * and the {@link ToolDispatcher} callback (coerced from a Groovy closure).
 */
class ToolDescriptorTest extends Specification {

    def 'should hold the descriptor fields'() {
        given:
        def inSchema = [type: 'object', properties: [name: [type: 'string']], required: ['name'], additionalProperties: false]
        def outSchema = [type: 'object', properties: [result: [type: 'string']], required: ['result'], additionalProperties: false]

        when:
        def desc = new ToolDescriptor('greet', 'Greet someone by name', inSchema, outSchema)

        then:
        desc.name == 'greet'
        desc.description == 'Greet someone by name'
        desc.inputSchema == inSchema
        desc.outputSchema == outSchema
    }

    def 'should coerce a closure to a ToolDispatcher and invoke it'() {
        given:
        String seenName = null
        String seenArgs = null
        ToolDispatcher dispatch = { String toolName, String argsJson ->
            seenName = toolName
            seenArgs = argsJson
            return '{"result":"ok"}'
        } as ToolDispatcher

        when:
        def out = dispatch.call('t', '{}')

        then:
        out == '{"result":"ok"}'
        seenName == 't'
        seenArgs == '{}'
    }
}
