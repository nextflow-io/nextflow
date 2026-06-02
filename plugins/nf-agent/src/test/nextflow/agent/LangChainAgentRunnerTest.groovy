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

import dev.langchain4j.data.message.AiMessage
import dev.langchain4j.data.message.ChatMessage
import dev.langchain4j.data.message.SystemMessage
import dev.langchain4j.data.message.UserMessage
import dev.langchain4j.model.chat.ChatModel
import dev.langchain4j.model.chat.request.json.JsonObjectSchema
import dev.langchain4j.model.chat.request.json.JsonSchema
import dev.langchain4j.model.chat.response.ChatResponse
import spock.lang.Specification

class LangChainAgentRunnerTest extends Specification {

    private static final Map ANSWER_SCHEMA = [
        type: 'object',
        properties: [answer: [type: 'string'], confidence: [type: 'number']],
        required: ['answer', 'confidence'],
        additionalProperties: false,
    ]

    def 'should compose prompt + input JSON, pass the schema, and return the assistant JSON'() {
        given:
        List<ChatMessage> captured = null
        // langchain4j ChatModel has no single abstract method (all default), so
        // mock it by overriding the chat(List<ChatMessage>) entry point used by the runner.
        ChatModel model = [
            chat: { List<ChatMessage> messages ->
                captured = messages
                ChatResponse.builder().aiMessage(AiMessage.from('{"answer":"ok","confidence":0.9}')).build()
            }
        ] as ChatModel

        and:
        JsonSchema capturedSchema = null
        def factory = Stub(ChatModelFactory) {
            createModel(_, _, _) >> { String id, int timeout, JsonSchema schema ->
                capturedSchema = schema
                model
            }
        }
        def runner = new LangChainAgentRunner(modelFactory: factory)
        def req = new AgentRunnerRequest(
            'openai/gpt-5-mini',
            'inst',
            'the prompt',
            5,
            [],
            ANSWER_SCHEMA,
            '{"text":"hi"}')

        when:
        def answer = runner.run(req)

        then:
        answer == '{"answer":"ok","confidence":0.9}'

        and: 'the schema was passed through to the model factory'
        capturedSchema != null
        capturedSchema.rootElement() instanceof JsonObjectSchema
        (capturedSchema.rootElement() as JsonObjectSchema).properties().containsKey('answer')

        and: 'the user message carries both the prompt and the input JSON'
        captured.size() == 2
        captured[0] instanceof SystemMessage
        (captured[0] as SystemMessage).text() == 'inst'
        captured[1] instanceof UserMessage
        def userText = (captured[1] as UserMessage).singleText()
        userText.contains('the prompt')
        userText.contains('{"text":"hi"}')
    }

    def 'should omit the system message and the input JSON when not provided'() {
        given:
        List<ChatMessage> captured = null
        ChatModel model = [
            chat: { List<ChatMessage> messages ->
                captured = messages
                ChatResponse.builder().aiMessage(AiMessage.from('{"answer":"ok"}')).build()
            }
        ] as ChatModel

        and:
        def runner = new LangChainAgentRunner(modelFactory: Stub(ChatModelFactory) {
            createModel(_, _, _) >> model
        })
        def req = new AgentRunnerRequest('openai/gpt-5-mini', null, 'just a prompt', 5, [], null, null)

        when:
        def answer = runner.run(req)

        then:
        answer == '{"answer":"ok"}'
        captured.size() == 1
        captured[0] instanceof UserMessage
        (captured[0] as UserMessage).singleText() == 'just a prompt'
    }

    def 'should fail when the model is missing'() {
        given:
        def runner = new LangChainAgentRunner()
        def req = new AgentRunnerRequest(null, 'inst', 'prompt', 5, [], null, null)

        when:
        runner.run(req)

        then:
        thrown(IllegalArgumentException)
    }

    def 'should pass the configured request timeout through to the model factory'() {
        given:
        ChatModel model = [
            chat: { List<ChatMessage> messages ->
                ChatResponse.builder().aiMessage(AiMessage.from('ok')).build()
            }
        ] as ChatModel

        and:
        int capturedTimeout = -1
        def factory = Stub(ChatModelFactory) {
            createModel(_, _, _) >> { String id, int timeout, JsonSchema schema ->
                capturedTimeout = timeout
                model
            }
        }
        def runner = new LangChainAgentRunner(modelFactory: factory)
        // 10th positional arg = requestTimeoutSeconds
        def req = new AgentRunnerRequest('openai/gpt-5-mini', null, 'prompt', 5, [], null, null, null, null, 90)

        when:
        runner.run(req)

        then:
        capturedTimeout == 90
    }

    def 'should fall back to the default request timeout when none is configured'() {
        given:
        ChatModel model = [
            chat: { List<ChatMessage> messages ->
                ChatResponse.builder().aiMessage(AiMessage.from('ok')).build()
            }
        ] as ChatModel

        and:
        int capturedTimeout = -1
        def factory = Stub(ChatModelFactory) {
            createModel(_, _, _) >> { String id, int timeout, JsonSchema schema ->
                capturedTimeout = timeout
                model
            }
        }
        def runner = new LangChainAgentRunner(modelFactory: factory)
        // requestTimeoutSeconds left at the default 0 -> built-in default (120)
        def req = new AgentRunnerRequest('openai/gpt-5-mini', null, 'prompt', 5, [], null, null)

        when:
        runner.run(req)

        then:
        capturedTimeout == 120
    }
}
