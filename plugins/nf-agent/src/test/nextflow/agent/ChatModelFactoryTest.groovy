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

import dev.langchain4j.model.chat.ChatModel
import dev.langchain4j.model.chat.request.json.JsonSchema
import spock.lang.Specification

class ChatModelFactoryTest extends Specification {

    def 'should split provider and model id'() {
        expect:
        ChatModelFactory.providerOf('openai/gpt-5-mini') == 'openai'
        ChatModelFactory.modelOf('openai/gpt-5-mini') == 'gpt-5-mini'
    }

    def 'should fail for a model id without a provider'() {
        when:
        ChatModelFactory.providerOf('gpt-5-mini')
        then:
        thrown(IllegalArgumentException)
    }

    def 'should fail for an unknown provider'() {
        when:
        new ChatModelFactory(apiKey: 'sk-test').createModel('acme/whatever', 30, null)
        then:
        def e = thrown(IllegalArgumentException)
        e.message.toLowerCase().contains('provider')
    }

    def 'should fail when the api key is missing'() {
        when:
        new ChatModelFactory(apiKey: null).createModel('openai/gpt-5-mini', 30, null)
        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains('OPENAI_API_KEY')
    }

    def 'should build an openai model when api key is present (no network call)'() {
        given:
        def factory = new ChatModelFactory(apiKey: 'sk-test')
        when:
        ChatModel model = factory.createModel('openai/gpt-5-mini', 30, null)
        then:
        model != null
    }

    def 'should build an openai model with a structured-output schema (no network call)'() {
        given:
        def factory = new ChatModelFactory(apiKey: 'sk-test')
        def schema = JsonSchemaMapper.toJsonSchema('Answer', [
            type: 'object',
            properties: [answer: [type: 'string']],
            required: ['answer'],
            additionalProperties: false,
        ])
        when:
        ChatModel model = factory.createModel('openai/gpt-5-mini', 30, schema)
        then:
        model != null
    }
}
