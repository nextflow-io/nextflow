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
import dev.langchain4j.data.message.UserMessage
import dev.langchain4j.model.chat.ChatModel
import dev.langchain4j.model.chat.request.ChatRequest
import dev.langchain4j.model.chat.request.json.JsonSchema
import dev.langchain4j.model.chat.response.ChatResponse
import dev.langchain4j.model.output.FinishReason
import spock.lang.Specification

/**
 * Exercises the M5 "final structuring turn": a tool-using (or skill-using) agent that
 * also declares a structured output. The tool loop runs schema-free (byte-for-byte the
 * M1-M4 behavior) and, only when {@code outputSchema != null}, a single stateless
 * structuring call converts the loop's free-text answer into schema-valid JSON.
 *
 * The {@link ChatModelFactory} is stubbed to hand back two different mock models
 * differentiated by the schema arg: a null schema returns the AiServices loop model
 * (overriding {@code chat(ChatRequest)}), a non-null schema returns the structuring
 * model (overriding {@code chat(List<ChatMessage>)}).
 */
class LangChainAgentToolStructuredTest extends Specification {

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

    /** A loop model that asks for the greet tool on the first turn, then answers free text. */
    private static ChatModel loopModel(Closure onCall = null) {
        int calls = 0
        return [
            chat: { ChatRequest req ->
                if( onCall != null ) onCall.call(req)
                calls++
                if( calls == 1 ) {
                    final ter = ToolExecutionRequest.builder()
                        .id('call-1').name('greet').arguments('{"name":"Ada"}').build()
                    return ChatResponse.builder().aiMessage(AiMessage.from([ter])).build()
                }
                return ChatResponse.builder().aiMessage(AiMessage.from('Hello Ada!')).build()
            }
        ] as ChatModel
    }

    def 'should run the tool loop schema-free then structure the final answer'() {
        given: 'a loop model (tool then free text) and a structuring model returning schema JSON'
        JsonSchema loopSchema = null
        JsonSchema structuringSchema = null
        List<ChatMessage> structuringMessages = null
        boolean structuringCalled = false

        def loop = loopModel()
        ChatModel structuring = [
            chat: { List<ChatMessage> messages ->
                structuringCalled = true
                structuringMessages = messages
                ChatResponse.builder().aiMessage(AiMessage.from('{"greeting":"Hello Ada!"}')).build()
            }
        ] as ChatModel

        and: 'a factory that routes on the schema arg (null -> loop, non-null -> structuring)'
        def factory = Stub(ChatModelFactory) {
            createModel(*_) >> { args ->
                final JsonSchema schema = args[2] as JsonSchema
                if( schema == null ) { loopSchema = schema; return loop }
                structuringSchema = schema
                return structuring
            }
        }
        def runner = new LangChainAgentRunner(modelFactory: factory)

        and: 'a stub dispatcher recording that the tool loop ran'
        List<List<String>> dispatched = []
        ToolDispatcher dispatch = { String name, String args ->
            dispatched << [name, args]; '{"greeting":"Hello Ada!"}'
        } as ToolDispatcher

        and: 'a tool request that ALSO carries a structured output schema'
        def descriptor = new ToolDescriptor('greet', 'greet someone', GREET_INPUT_SCHEMA, GREET_OUTPUT_SCHEMA)
        def req = new AgentRunnerRequest(
            model: 'openai/gpt-5-mini', instruction: 'Use the greet tool.', prompt: 'greet Ada',
            maxIterations: 5, tools: [], outputSchema: GREET_OUTPUT_SCHEMA, inputJson: null,
            toolSpecs: [descriptor], dispatch: dispatch, requestTimeoutSeconds: 0, agentName: 'assistant')

        when:
        def answer = runner.run(req)

        then: 'the tool loop actually ran'
        dispatched == [['greet', '{"name":"Ada"}']]

        and: 'the loop model received a NULL schema (schema-free loop preserved)'
        loopSchema == null

        and: 'the structuring model received the non-null schema'
        structuringCalled
        structuringSchema != null

        and: 'the returned value is the structuring JSON, not the free text'
        answer == '{"greeting":"Hello Ada!"}'

        and: 'the structuring conversation is exactly a system instruction + the loop answer as user text'
        structuringMessages.size() == 2
        structuringMessages[0] instanceof SystemMessage
        (structuringMessages[0] as SystemMessage).text().contains('structured JSON')
        structuringMessages[1] instanceof UserMessage
        (structuringMessages[1] as UserMessage).singleText() == 'Hello Ada!'
    }

    def 'should not make a structuring turn when no outputSchema is set'() {
        given: 'a loop model and a structuring model that must never be called'
        boolean structuringCalled = false
        def loop = loopModel()
        ChatModel structuring = [
            chat: { List<ChatMessage> messages ->
                structuringCalled = true
                ChatResponse.builder().aiMessage(AiMessage.from('{}')).build()
            }
        ] as ChatModel
        def factory = Stub(ChatModelFactory) {
            createModel(*_) >> { args -> (args[2] as JsonSchema) == null ? loop : structuring }
        }
        def runner = new LangChainAgentRunner(modelFactory: factory)
        ToolDispatcher dispatch = { String n, String a -> '{"greeting":"Hello Ada!"}' } as ToolDispatcher
        def descriptor = new ToolDescriptor('greet', 'greet someone', GREET_INPUT_SCHEMA, null)
        def req = new AgentRunnerRequest(
            model: 'openai/gpt-5-mini', instruction: 'Use the greet tool.', prompt: 'greet Ada',
            maxIterations: 5, tools: [], outputSchema: null, inputJson: null,
            toolSpecs: [descriptor], dispatch: dispatch, requestTimeoutSeconds: 0, agentName: 'assistant')

        when:
        def answer = runner.run(req)

        then: 'the free-text answer is returned verbatim and no structuring turn happened'
        answer == 'Hello Ada!'
        !structuringCalled
    }

    def 'should structure a skills-only structured request'() {
        given: 'a skills-only request (no toolSpecs) that answers immediately, plus a structuring model'
        def loop = [
            chat: { ChatRequest req ->
                ChatResponse.builder().aiMessage(AiMessage.from('Hi Ada')).build()
            }
        ] as ChatModel
        boolean structuringCalled = false
        ChatModel structuring = [
            chat: { List<ChatMessage> messages ->
                structuringCalled = true
                ChatResponse.builder().aiMessage(AiMessage.from('{"greeting":"Hi Ada"}')).build()
            }
        ] as ChatModel
        def factory = Stub(ChatModelFactory) {
            createModel(*_) >> { args -> (args[2] as JsonSchema) == null ? loop : structuring }
        }
        def runner = new LangChainAgentRunner(modelFactory: factory)
        def skill = new SkillDescriptor('greet', 'a greeting skill', 'say hi politely', [])
        def req = new AgentRunnerRequest(
            model: 'openai/gpt-5-mini', instruction: 'You are helpful.', prompt: 'greet Ada',
            maxIterations: 5, tools: [], outputSchema: GREET_OUTPUT_SCHEMA, inputJson: null,
            toolSpecs: null, dispatch: null, requestTimeoutSeconds: 0, goal: null,
            agentName: 'assistant', skills: [skill])

        when:
        def answer = runner.run(req)

        then: 'the skills path (runWithTools) also runs the structuring turn'
        structuringCalled
        answer == '{"greeting":"Hi Ada"}'
    }

    def 'should throw a refusal-flavored error when the structuring turn is content-filtered'() {
        given:
        def loop = loopModel()
        ChatModel structuring = [
            chat: { List<ChatMessage> messages ->
                ChatResponse.builder().aiMessage(AiMessage.from('')).finishReason(FinishReason.CONTENT_FILTER).build()
            }
        ] as ChatModel
        def factory = Stub(ChatModelFactory) {
            createModel(*_) >> { args -> (args[2] as JsonSchema) == null ? loop : structuring }
        }
        def runner = new LangChainAgentRunner(modelFactory: factory)
        ToolDispatcher dispatch = { String n, String a -> '{"greeting":"Hello Ada!"}' } as ToolDispatcher
        def descriptor = new ToolDescriptor('greet', 'greet someone', GREET_INPUT_SCHEMA, GREET_OUTPUT_SCHEMA)
        def req = new AgentRunnerRequest(
            model: 'openai/gpt-5-mini', instruction: 'Use the greet tool.', prompt: 'greet Ada',
            maxIterations: 5, tools: [], outputSchema: GREET_OUTPUT_SCHEMA, inputJson: null,
            toolSpecs: [descriptor], dispatch: dispatch, requestTimeoutSeconds: 0, agentName: 'assistant')

        when:
        runner.run(req)

        then:
        def e = thrown(AgentStructuredOutputException)
        e.message.contains('assistant')
    }

    def 'should throw a refusal-flavored error when the structuring turn returns blank content'() {
        given:
        def loop = loopModel()
        ChatModel structuring = [
            chat: { List<ChatMessage> messages ->
                // blank text, normal STOP finish reason -> must NOT fall through to a blank result
                ChatResponse.builder().aiMessage(AiMessage.from('   ')).finishReason(FinishReason.STOP).build()
            }
        ] as ChatModel
        def factory = Stub(ChatModelFactory) {
            createModel(*_) >> { args -> (args[2] as JsonSchema) == null ? loop : structuring }
        }
        def runner = new LangChainAgentRunner(modelFactory: factory)
        ToolDispatcher dispatch = { String n, String a -> '{"greeting":"Hello Ada!"}' } as ToolDispatcher
        def descriptor = new ToolDescriptor('greet', 'greet someone', GREET_INPUT_SCHEMA, GREET_OUTPUT_SCHEMA)
        def req = new AgentRunnerRequest(
            model: 'openai/gpt-5-mini', instruction: 'Use the greet tool.', prompt: 'greet Ada',
            maxIterations: 5, tools: [], outputSchema: GREET_OUTPUT_SCHEMA, inputJson: null,
            toolSpecs: [descriptor], dispatch: dispatch, requestTimeoutSeconds: 0, agentName: 'assistant')

        when:
        runner.run(req)

        then:
        thrown(AgentStructuredOutputException)
    }

    def 'should thread the request timeout into the structuring model'() {
        given: 'capture the timeout passed for the structuring (non-null schema) call'
        int structuringTimeout = -1
        def loop = loopModel()
        ChatModel structuring = [
            chat: { List<ChatMessage> messages ->
                ChatResponse.builder().aiMessage(AiMessage.from('{"greeting":"Hello Ada!"}')).build()
            }
        ] as ChatModel
        def factory = Stub(ChatModelFactory) {
            createModel(*_) >> { args ->
                final JsonSchema schema = args[2] as JsonSchema
                if( schema != null ) structuringTimeout = args[1] as int
                schema == null ? loop : structuring
            }
        }
        def runner = new LangChainAgentRunner(modelFactory: factory)
        ToolDispatcher dispatch = { String n, String a -> '{"greeting":"Hello Ada!"}' } as ToolDispatcher
        def descriptor = new ToolDescriptor('greet', 'greet someone', GREET_INPUT_SCHEMA, GREET_OUTPUT_SCHEMA)
        def req = new AgentRunnerRequest(
            model: 'openai/gpt-5-mini', instruction: 'Use the greet tool.', prompt: 'greet Ada',
            maxIterations: 5, tools: [], outputSchema: GREET_OUTPUT_SCHEMA, inputJson: null,
            toolSpecs: [descriptor], dispatch: dispatch, requestTimeoutSeconds: 90, agentName: 'assistant')

        when:
        runner.run(req)

        then: 'the per-request timeout was applied to the structuring call'
        structuringTimeout == 90
    }
}
