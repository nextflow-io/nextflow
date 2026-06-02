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

import dev.langchain4j.agent.tool.ToolExecutionRequest
import dev.langchain4j.data.message.AiMessage
import dev.langchain4j.data.message.ChatMessage
import dev.langchain4j.data.message.ToolExecutionResultMessage
import dev.langchain4j.data.message.UserMessage
import dev.langchain4j.model.chat.ChatModel
import dev.langchain4j.model.chat.request.ChatRequest
import dev.langchain4j.model.chat.request.json.JsonSchema
import dev.langchain4j.model.chat.response.ChatResponse
import spock.lang.Specification

class LangChainAgentToolLoopTest extends Specification {

    private static final Map GREET_INPUT_SCHEMA = [
        type: 'object',
        properties: [name: [type: 'string']],
        required: ['name'],
        additionalProperties: false,
    ]

    private static final Map GREET_OUTPUT_SCHEMA = [
        type: 'object',
        properties: [greeting: [type: 'string']],
        required: ['greeting'],
        additionalProperties: false,
    ]

    def 'should drive the tool loop: call dispatch then return the final text'() {
        given: 'a mock model that requests the greet tool first, then answers'
        List<ChatRequest> capturedRequests = []
        int calls = 0
        // ChatModel has no single abstract method (all chat overloads are default);
        // the tool loop calls chat(ChatRequest), so override that overload.
        ChatModel model = [
            chat: { ChatRequest req ->
                capturedRequests << req
                calls++
                if( calls == 1 ) {
                    // first turn: ask to run the greet tool
                    final ter = ToolExecutionRequest.builder()
                        .id('call-1')
                        .name('greet')
                        .arguments('{"name":"Ada"}')
                        .build()
                    return ChatResponse.builder().aiMessage(AiMessage.from([ter])).build()
                }
                // second turn: final plain-text answer (no tool requests)
                return ChatResponse.builder().aiMessage(AiMessage.from('Final: greeted Ada')).build()
            }
        ] as ChatModel

        and: 'a stub dispatcher recording the call and returning a canned JSON result'
        List<List<String>> dispatched = []
        ToolDispatcher dispatch = { String name, String args ->
            dispatched << [name, args]
            return '{"greeting":"Hello Ada!"}'
        } as ToolDispatcher

        and: 'the runner wired to a factory that never forces a responseFormat'
        JsonSchema capturedSchema = null
        boolean factoryCalled = false
        def factory = Stub(ChatModelFactory) {
            createModel(_, _, _) >> { String id, int timeout, JsonSchema schema ->
                factoryCalled = true
                capturedSchema = schema
                model
            }
        }
        def runner = new LangChainAgentRunner(modelFactory: factory)

        and: 'a request carrying the greet tool spec and the dispatch callback'
        def descriptor = new ToolDescriptor('greet', 'greet someone', GREET_INPUT_SCHEMA, GREET_OUTPUT_SCHEMA)
        def req = new AgentRunnerRequest(
            'openai/gpt-5-mini',
            'Use the greet tool.',
            'greet Ada',
            5,
            [],
            null,
            null,
            [descriptor],
            dispatch)

        when:
        def answer = runner.run(req)

        then: 'the final text is returned'
        answer == 'Final: greeted Ada'

        and: 'the dispatcher was called exactly once with the right name and args'
        dispatched == [['greet', '{"name":"Ada"}']]

        and: 'the model was chatted twice'
        calls == 2

        and: 'tools were advertised on every request, including the greet spec'
        capturedRequests.size() == 2
        capturedRequests.every { it.toolSpecifications() && it.toolSpecifications()*.name().contains('greet') }

        and: 'no structured-output schema was forced (tools XOR responseFormat)'
        factoryCalled
        capturedSchema == null

        and: 'the second request carried the tool result message after the assistant turn'
        List<ChatMessage> second = capturedRequests[1].messages()
        def resultMsg = second.find { it instanceof ToolExecutionResultMessage } as ToolExecutionResultMessage
        resultMsg != null
        resultMsg.toolName() == 'greet'
        resultMsg.text() == '{"greeting":"Hello Ada!"}'

        and: 'the second request still carried the original user prompt'
        second.find { it instanceof UserMessage && (it as UserMessage).singleText().contains('greet Ada') } != null
    }

    def 'should throw when the iteration cap is exceeded without a final answer'() {
        given: 'a model that always requests the tool, never answering'
        ChatModel model = [
            chat: { ChatRequest req ->
                final ter = ToolExecutionRequest.builder()
                    .id('loop')
                    .name('greet')
                    .arguments('{"name":"Ada"}')
                    .build()
                ChatResponse.builder().aiMessage(AiMessage.from([ter])).build()
            }
        ] as ChatModel

        and:
        int dispatchCalls = 0
        ToolDispatcher dispatch = { String name, String args -> dispatchCalls++; '{"greeting":"Hello"}' } as ToolDispatcher
        def runner = new LangChainAgentRunner(modelFactory: Stub(ChatModelFactory) {
            createModel(_, _, _) >> model
        })
        def descriptor = new ToolDescriptor('greet', 'greet someone', GREET_INPUT_SCHEMA, null)
        def req = new AgentRunnerRequest('openai/gpt-5-mini', null, 'greet Ada', 2, [], null, null, [descriptor], dispatch)

        when:
        runner.run(req)

        then:
        thrown(IllegalStateException)

        and: 'the dispatcher ran once per iteration up to the cap'
        dispatchCalls == 2
    }
}
