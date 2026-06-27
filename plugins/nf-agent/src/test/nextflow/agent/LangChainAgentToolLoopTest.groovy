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
import dev.langchain4j.data.message.SystemMessage
import dev.langchain4j.data.message.ToolExecutionResultMessage
import dev.langchain4j.data.message.UserMessage
import dev.langchain4j.model.chat.ChatModel
import dev.langchain4j.model.chat.request.ChatRequest
import dev.langchain4j.model.chat.request.json.JsonSchema
import dev.langchain4j.model.chat.response.ChatResponse
import spock.lang.Specification

/**
 * Exercises {@link LangChainAgentRunner#runWithTools} which is now driven by a
 * langchain4j {@code AiServices} proxy. AiServices routes the LLM call through
 * {@code ChatModel.chat(ChatRequest)}, so a Groovy-closure-coerced ChatModel
 * mock overriding that overload remains the injection seam. AiServices builds
 * each ChatRequest from the shared (seeded) chat memory, so each request's
 * {@code messages()} snapshot reflects the accumulating memory contents.
 */
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

    def 'should drive the AiServices tool loop: call dispatch then return the final text'() {
        given: 'a mock model that requests the greet tool first, then answers'
        List<ChatRequest> capturedRequests = []
        int calls = 0
        // AiServices invokes chat(ChatRequest); override only that overload.
        ChatModel model = [
            chat: { ChatRequest req ->
                // snapshot the memory-derived messages as the proxy sees them
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

        and: 'the model was chatted twice (one round-trip per tool turn plus the final)'
        calls == 2

        and: 'tools were advertised on every request, including the greet spec'
        capturedRequests.size() == 2
        capturedRequests.every { it.toolSpecifications() && it.toolSpecifications()*.name().contains('greet') }

        and: 'no structured-output schema was forced (tools XOR responseFormat)'
        factoryCalled
        capturedSchema == null

        and: 'the FIRST request memory held exactly: system message + one user message (full prompt)'
        List<ChatMessage> first = capturedRequests[0].messages()
        first.count { it instanceof SystemMessage } == 1
        first.count { it instanceof UserMessage } == 1
        (first.find { it instanceof SystemMessage } as SystemMessage).text() == 'Use the greet tool.'
        (first.find { it instanceof UserMessage } as UserMessage).singleText() == 'greet Ada'

        and: 'the SECOND request memory grew with the assistant tool-request turn and the tool result'
        List<ChatMessage> second = capturedRequests[1].messages()
        // still exactly one user message — no prompt duplication
        second.count { it instanceof UserMessage } == 1
        // the assistant turn carrying the tool request is present
        second.find { it instanceof AiMessage && (it as AiMessage).hasToolExecutionRequests() } != null
        // the tool result message was fed back with the right name and JSON payload
        def resultMsg = second.find { it instanceof ToolExecutionResultMessage } as ToolExecutionResultMessage
        resultMsg != null
        resultMsg.toolName() == 'greet'
        resultMsg.text() == '{"greeting":"Hello Ada!"}'

        and: 'the original user prompt persisted unchanged across the loop'
        (second.find { it instanceof UserMessage } as UserMessage).singleText() == 'greet Ada'
    }

    def 'should compose user text with the input JSON as a single user message'() {
        given: 'a model that answers immediately on the first turn'
        List<ChatRequest> capturedRequests = []
        ChatModel model = [
            chat: { ChatRequest req ->
                capturedRequests << req
                return ChatResponse.builder().aiMessage(AiMessage.from('done')).build()
            }
        ] as ChatModel

        and:
        ToolDispatcher dispatch = { String name, String args -> '{}' } as ToolDispatcher
        def factory = Stub(ChatModelFactory) { createModel(_, _, _) >> model }
        def runner = new LangChainAgentRunner(modelFactory: factory)
        def descriptor = new ToolDescriptor('greet', 'greet someone', GREET_INPUT_SCHEMA, null)
        def req = new AgentRunnerRequest(
            'openai/gpt-5-mini',
            'sys',
            'do it',
            5,
            [],
            null,
            '{"k":"v"}',
            [descriptor],
            dispatch)

        when:
        def answer = runner.run(req)

        then: 'the answer is returned'
        answer == 'done'

        and: 'exactly one user message carrying prompt + input JSON'
        List<ChatMessage> msgs = capturedRequests[0].messages()
        msgs.count { it instanceof UserMessage } == 1
        (msgs.find { it instanceof UserMessage } as UserMessage).singleText() == 'do it\n\nInput (JSON):\n{"k":"v"}'
    }

    def 'should succeed with exactly one tool round-trip when maxIterations=1 (verifies the +1 cap offset)'() {
        // Rationale: maxIterations=1 → maxSequentialToolsInvocations is passed as 2.
        // A cap of 2 permits one tool round-trip (tool call + final answer).
        // If the code wrongly passed maxIterations (=1) the cap would be 1, which
        // allows zero round-trips, and the run would throw before the tool executes.
        // This test therefore fails if the "+1" in maxIterations+1 is dropped.
        given: 'a model that requests the greet tool on the first call, then answers'
        int calls = 0
        ChatModel model = [
            chat: { ChatRequest req ->
                calls++
                if( calls == 1 ) {
                    final ter = ToolExecutionRequest.builder()
                        .id('call-boundary')
                        .name('greet')
                        .arguments('{"name":"Boundary"}')
                        .build()
                    return ChatResponse.builder().aiMessage(AiMessage.from([ter])).build()
                }
                // second turn: final plain-text answer
                return ChatResponse.builder().aiMessage(AiMessage.from('Boundary answer')).build()
            }
        ] as ChatModel

        and:
        int dispatchCalls = 0
        ToolDispatcher dispatch = { String name, String args -> dispatchCalls++; '{"greeting":"Hello Boundary!"}' } as ToolDispatcher
        def runner = new LangChainAgentRunner(modelFactory: Stub(ChatModelFactory) {
            createModel(_, _, _) >> model
        })
        def descriptor = new ToolDescriptor('greet', 'greet someone', GREET_INPUT_SCHEMA, null)
        // maxIterations=1: only one tool round-trip is permitted
        def req = new AgentRunnerRequest('openai/gpt-5-mini', null, 'greet Boundary', 1, [], null, null, [descriptor], dispatch)

        when:
        def answer = runner.run(req)

        then: 'the run succeeds and returns the final text'
        answer == 'Boundary answer'

        and: 'the dispatcher was called exactly once'
        dispatchCalls == 1
    }

    def 'should throw IllegalStateException when the iteration cap is exceeded without a final answer'() {
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

        then: 'the cap RuntimeException is surfaced as the historical IllegalStateException'
        def e = thrown(IllegalStateException)
        e.message == 'Agent exceeded the maximum number of tool-call iterations (2)'
    }

    def 'should fold goal into the single system message (tool path)'() {
        given:
        List<ChatRequest> capturedRequests = []
        ChatModel model = [
            chat: { ChatRequest req ->
                capturedRequests << req
                ChatResponse.builder()
                    .aiMessage(AiMessage.from('done')).build()
            }
        ] as ChatModel
        ToolDispatcher dispatch = { String n, String a -> '{}' } as ToolDispatcher
        def factory = Stub(ChatModelFactory) { createModel(_, _, _) >> model }
        def runner = new LangChainAgentRunner(modelFactory: factory)
        def descriptor = new ToolDescriptor('greet', 'greet', [type:'object', properties:[name:[type:'string']], required:['name'], additionalProperties:false], null)
        def req = new AgentRunnerRequest(
            model: 'openai/gpt-5-mini', instruction: 'You are careful.', prompt: 'go',
            maxIterations: 5, tools: [], outputSchema: null, inputJson: null,
            toolSpecs: [descriptor], dispatch: dispatch, requestTimeoutSeconds: 0,
            goal: 'assemble then QC')

        when:
        def answer = runner.run(req)

        then:
        answer == 'done'
        def sys = capturedRequests[0].messages().find { it instanceof SystemMessage } as SystemMessage
        capturedRequests[0].messages().count { it instanceof SystemMessage } == 1
        sys.text().contains('You are careful.')
        sys.text().contains('assemble then QC')
    }

    def 'should produce a system message from goal alone (no instruction)'() {
        given:
        List<ChatRequest> capturedRequests = []
        ChatModel model = [
            chat: { ChatRequest req ->
                capturedRequests << req
                ChatResponse.builder()
                    .aiMessage(AiMessage.from('ok')).build()
            }
        ] as ChatModel
        ToolDispatcher dispatch = { String n, String a -> '{}' } as ToolDispatcher
        def runner = new LangChainAgentRunner(modelFactory: Stub(ChatModelFactory) { createModel(_, _, _) >> model })
        def descriptor = new ToolDescriptor('greet', 'greet', [type:'object', properties:[name:[type:'string']], required:['name'], additionalProperties:false], null)
        def req = new AgentRunnerRequest(
            model: 'openai/gpt-5-mini', instruction: null, prompt: 'go',
            maxIterations: 5, tools: [], outputSchema: null, inputJson: null,
            toolSpecs: [descriptor], dispatch: dispatch, requestTimeoutSeconds: 0,
            goal: 'reach the objective')

        when:
        runner.run(req)

        then:
        def msgs = capturedRequests[0].messages()
        msgs.count { it instanceof SystemMessage } == 1
        (msgs.find { it instanceof SystemMessage } as SystemMessage).text().contains('reach the objective')
    }

    def 'should seed no system message when neither instruction nor goal is set'() {
        given:
        List<ChatRequest> capturedRequests = []
        ChatModel model = [
            chat: { ChatRequest req ->
                capturedRequests << req
                ChatResponse.builder()
                    .aiMessage(AiMessage.from('x')).build()
            }
        ] as ChatModel
        ToolDispatcher dispatch = { String n, String a -> '{}' } as ToolDispatcher
        def runner = new LangChainAgentRunner(modelFactory: Stub(ChatModelFactory) { createModel(_, _, _) >> model })
        def descriptor = new ToolDescriptor('greet', 'greet', [type:'object', properties:[name:[type:'string']], required:['name'], additionalProperties:false], null)
        def req = new AgentRunnerRequest(
            model: 'openai/gpt-5-mini', instruction: null, prompt: 'go',
            maxIterations: 5, tools: [], outputSchema: null, inputJson: null,
            toolSpecs: [descriptor], dispatch: dispatch, requestTimeoutSeconds: 0, goal: null)

        when:
        runner.run(req)

        then:
        capturedRequests[0].messages().count { it instanceof SystemMessage } == 0
    }

    def 'should add the reasoning-narration directive to the system message when tracing'() {
        given:
        List<ChatRequest> capturedRequests = []
        ChatModel model = [
            chat: { ChatRequest req ->
                capturedRequests << req
                ChatResponse.builder().aiMessage(AiMessage.from('done')).build()
            }
        ] as ChatModel
        ToolDispatcher dispatch = { String n, String a -> '{}' } as ToolDispatcher
        // tracing builds the model via the 4-arg createModel (with listeners)
        def runner = new LangChainAgentRunner(modelFactory: Stub(ChatModelFactory) { createModel(_, _, _, _) >> model })
        def descriptor = new ToolDescriptor('greet', 'greet', [type:'object', properties:[name:[type:'string']], required:['name'], additionalProperties:false], null)
        def req = new AgentRunnerRequest(
            model: 'openai/gpt-5-mini', instruction: 'You are careful.', prompt: 'go',
            maxIterations: 5, tools: [], outputSchema: null, inputJson: null,
            toolSpecs: [descriptor], dispatch: dispatch, requestTimeoutSeconds: 0,
            goal: null, agentName: 't', trace: true)

        when:
        runner.run(req)

        then:
        def sys = capturedRequests[0].messages().find { it instanceof SystemMessage } as SystemMessage
        sys.text().contains('You are careful.')
        sys.text().contains('briefly state your reasoning')
    }

    def 'should not add the narration directive when not tracing'() {
        given:
        List<ChatRequest> capturedRequests = []
        ChatModel model = [
            chat: { ChatRequest req ->
                capturedRequests << req
                ChatResponse.builder().aiMessage(AiMessage.from('done')).build()
            }
        ] as ChatModel
        ToolDispatcher dispatch = { String n, String a -> '{}' } as ToolDispatcher
        def runner = new LangChainAgentRunner(modelFactory: Stub(ChatModelFactory) { createModel(_, _, _) >> model })
        def descriptor = new ToolDescriptor('greet', 'greet', [type:'object', properties:[name:[type:'string']], required:['name'], additionalProperties:false], null)
        def req = new AgentRunnerRequest(
            model: 'openai/gpt-5-mini', instruction: 'You are careful.', prompt: 'go',
            maxIterations: 5, tools: [], outputSchema: null, inputJson: null,
            toolSpecs: [descriptor], dispatch: dispatch, requestTimeoutSeconds: 0,
            goal: null)

        when:
        runner.run(req)

        then:
        def sys = capturedRequests[0].messages().find { it instanceof SystemMessage } as SystemMessage
        sys.text().contains('You are careful.')
        !sys.text().contains('briefly state your reasoning')
    }
}
