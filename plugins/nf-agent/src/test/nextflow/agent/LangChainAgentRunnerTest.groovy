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
import dev.langchain4j.model.chat.response.ChatResponseMetadata
import spock.lang.Specification

class LangChainAgentRunnerTest extends Specification {

    def cleanup() {
        AgentCallInfo.clear()
    }

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
            // createModel(modelId, timeout, schema, temperature, apiKey, baseUrl[, listeners])
            createModel(*_) >> { args ->
                capturedSchema = args[2] as JsonSchema
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
            createModel(*_) >> model
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

    def 'should fold goal into the system message on the single-shot path'() {
        given:
        List<ChatMessage> captured = null
        ChatModel model = [
            chat: { List<ChatMessage> messages ->
                captured = messages
                ChatResponse.builder().aiMessage(AiMessage.from('result')).build()
            }
        ] as ChatModel

        and:
        def runner = new LangChainAgentRunner(modelFactory: Stub(ChatModelFactory) {
            createModel(*_) >> model
        })
        def req = new AgentRunnerRequest(
            model: 'openai/gpt-5-mini', instruction: 'Be helpful.', prompt: 'do it',
            maxIterations: 5, tools: [], outputSchema: null, inputJson: null,
            toolSpecs: null, dispatch: null, requestTimeoutSeconds: 0,
            goal: 'reach the summit')

        when:
        def answer = runner.run(req)

        then:
        answer == 'result'

        and: 'exactly one system message carrying both instruction and goal'
        captured.count { it instanceof SystemMessage } == 1
        def sys = captured.find { it instanceof SystemMessage } as SystemMessage
        sys.text().contains('Be helpful.')
        sys.text().contains('reach the summit')

        and: 'the user message follows'
        captured.count { it instanceof UserMessage } == 1
        (captured.find { it instanceof UserMessage } as UserMessage).singleText() == 'do it'
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
            createModel(*_) >> { args ->
                capturedTimeout = args[1] as int
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
            createModel(*_) >> { args ->
                capturedTimeout = args[1] as int
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

    def 'should thread the request temperature through on the single-shot path'() {
        given:
        ChatModel model = [
            chat: { List<ChatMessage> messages ->
                ChatResponse.builder().aiMessage(AiMessage.from('ok')).build()
            }
        ] as ChatModel

        and:
        Double capturedTemp = -999d
        def factory = Stub(ChatModelFactory) {
            createModel(*_) >> { args ->
                capturedTemp = args[3] as Double
                model
            }
        }
        def runner = new LangChainAgentRunner(modelFactory: factory)
        // a tool-free single-shot request carrying a pinned temperature (task path pins 0.0)
        def req = new AgentRunnerRequest(model: 'openai/gpt-5-mini', prompt: 'prompt', temperature: 0.0d)

        when:
        runner.run(req)

        then: 'the request temperature is threaded through (boxed, non-null) to the model factory'
        capturedTemp == 0.0d
        capturedTemp instanceof Double
    }

    def 'should capture the resolved model snapshot into AgentCallInfo after model.chat (design §9.5/D6)'() {
        given: 'a model whose response metadata carries the concrete snapshot'
        ChatModel model = [
            chat: { List<ChatMessage> messages ->
                ChatResponse.builder()
                    .aiMessage(AiMessage.from('ok'))
                    .metadata(ChatResponseMetadata.builder().modelName('gpt-4o-2024-08-06').build())
                    .build()
            }
        ] as ChatModel
        def runner = new LangChainAgentRunner(modelFactory: Stub(ChatModelFactory) {
            createModel(*_) >> model
        })
        def req = new AgentRunnerRequest(model: 'openai/gpt-4o', prompt: 'p')

        when:
        runner.run(req)

        then: 'the concrete snapshot crossed the SPI via the core-owned ThreadLocal'
        AgentCallInfo.consumeResolvedModel() == 'gpt-4o-2024-08-06'
    }

    def 'should leave temperature unset (null) on the legacy tools path'() {
        given: 'a request advertising a tool so the runner takes the runWithTools path'
        Double capturedTemp = -999d
        ChatModel model = [
            chat: { List<ChatMessage> messages ->
                ChatResponse.builder().aiMessage(AiMessage.from('ok')).build()
            }
        ] as ChatModel
        def factory = Stub(ChatModelFactory) {
            createModel(*_) >> { args ->
                // buildChatModel passes (modelId, timeout, null, (Double)null[, listeners])
                capturedTemp = args[3] as Double
                model
            }
        }
        def runner = new LangChainAgentRunner(modelFactory: factory)
        def descriptor = new ToolDescriptor('noop', 'a tool', [type: 'object', properties: [:], required: [], additionalProperties: false], [:])
        def dispatch = { String name, String argsJson -> '{}' } as ToolDispatcher
        def req = new AgentRunnerRequest(
            model: 'openai/gpt-5-mini', prompt: 'prompt', maxIterations: 3,
            toolSpecs: [descriptor], dispatch: dispatch)

        when: 'the tool loop is driven (the AiServices call may fail on the bare model stub - irrelevant here)'
        try { runner.run(req) } catch( Throwable ignored ) {}

        then: 'buildChatModel invoked createModel with a null temperature (legacy behavior unchanged)'
        capturedTemp == null
    }

    def 'should thread the resolved credential and endpoint through on the single-shot path'() {
        given:
        ChatModel model = [
            chat: { List<ChatMessage> messages ->
                ChatResponse.builder().aiMessage(AiMessage.from('ok')).build()
            }
        ] as ChatModel

        and: 'the runner must not read the environment: both values come from the request'
        String capturedKey = null
        String capturedUrl = null
        def factory = Stub(ChatModelFactory) {
            createModel(*_) >> { args ->
                capturedKey = args[4] as String
                capturedUrl = args[5] as String
                model
            }
        }
        def runner = new LangChainAgentRunner(modelFactory: factory)
        def req = new AgentRunnerRequest(
            model: 'openai/gpt-5-mini', prompt: 'prompt',
            apiKey: 'sk-from-config', baseUrl: 'http://localhost:8000/v1')

        when:
        runner.run(req)

        then:
        capturedKey == 'sk-from-config'
        capturedUrl == 'http://localhost:8000/v1'
    }

    def 'should thread the resolved credential and endpoint through on the tools path'() {
        given: 'a request advertising a tool so the runner takes the runWithTools path'
        String capturedKey = null
        String capturedUrl = null
        ChatModel model = [
            chat: { List<ChatMessage> messages ->
                ChatResponse.builder().aiMessage(AiMessage.from('ok')).build()
            }
        ] as ChatModel
        def factory = Stub(ChatModelFactory) {
            createModel(*_) >> { args ->
                capturedKey = args[4] as String
                capturedUrl = args[5] as String
                model
            }
        }
        def runner = new LangChainAgentRunner(modelFactory: factory)
        def descriptor = new ToolDescriptor('noop', 'a tool', [type: 'object', properties: [:], required: [], additionalProperties: false], [:])
        def dispatch = { String name, String argsJson -> '{}' } as ToolDispatcher
        def req = new AgentRunnerRequest(
            model: 'openai/gpt-5-mini', prompt: 'prompt', maxIterations: 3,
            toolSpecs: [descriptor], dispatch: dispatch,
            apiKey: 'sk-from-config', baseUrl: 'http://localhost:8000/v1')

        when: 'the tool loop is driven (the AiServices call may fail on the bare model stub - irrelevant here)'
        try { runner.run(req) } catch( Throwable ignored ) {}

        then: 'buildChatModel forwarded both, before the trailing listeners argument'
        capturedKey == 'sk-from-config'
        capturedUrl == 'http://localhost:8000/v1'
    }
}
