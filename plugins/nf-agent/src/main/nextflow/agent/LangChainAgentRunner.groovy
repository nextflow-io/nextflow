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
import dev.langchain4j.agent.tool.ToolSpecification
import dev.langchain4j.data.message.AiMessage
import dev.langchain4j.data.message.ChatMessage
import dev.langchain4j.data.message.SystemMessage
import dev.langchain4j.data.message.ToolExecutionResultMessage
import dev.langchain4j.data.message.UserMessage
import dev.langchain4j.model.chat.ChatModel
import dev.langchain4j.model.chat.request.ChatRequest
import dev.langchain4j.model.chat.request.json.JsonSchema
import dev.langchain4j.model.chat.response.ChatResponse
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import org.pf4j.Extension

/**
 * langchain4j-backed {@link AgentRunner}.
 *
 * When the request declares no tools, the runner does a single-shot chat
 * (optionally constraining the model to a structured-output JSON schema). When
 * the request declares one or more tools, the runner instead drives a manual
 * tool-call loop: it advertises the tool specifications on each chat request and,
 * for every tool-execution request the model emits, dispatches the call back into
 * core (which runs the real module) and feeds the result to the model, looping
 * until the model returns a final text answer or the iteration cap is reached.
 *
 * For Phase 2 tools and structured-output are mutually exclusive: when tools are
 * in play no {@code responseFormat} is forced on the model.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@Extension
@CompileStatic
class LangChainAgentRunner implements AgentRunner {

    private static final int DEFAULT_TIMEOUT_SECONDS = 120

    private static final int DEFAULT_MAX_ITERATIONS = 10

    ChatModelFactory modelFactory = new ChatModelFactory()

    @Override
    String run(AgentRunnerRequest request) {
        if( !request.model )
            throw new IllegalArgumentException("Agent `model` directive is required")

        return request.toolSpecs
            ? runWithTools(request)
            : runSingleShot(request)
    }

    /**
     * Single-shot chat (no tool calls). When an output record type is declared
     * the model is constrained to a structured-output JSON schema.
     */
    private String runSingleShot(AgentRunnerRequest request) {
        // derive the structured-output schema (when an output record type is declared)
        final JsonSchema schema = request.outputSchema
            ? JsonSchemaMapper.toJsonSchema('Output', request.outputSchema)
            : null

        final model = modelFactory.createModel(request.model, DEFAULT_TIMEOUT_SECONDS, schema)

        final List<ChatMessage> messages = composeMessages(request)

        log.debug "Running agent model=${request.model}; messages=${messages.size()}; structured=${schema != null}"
        final ChatResponse response = model.chat(messages)
        return response.aiMessage().text()
    }

    /**
     * Manual tool-call loop. The tool specifications are advertised on every
     * chat request; for each tool-execution request the model emits, the dispatch
     * callback runs the real module and its result is fed back to the model. The
     * loop ends when the model returns a final text answer (no tool requests) or
     * the iteration cap is reached.
     */
    private String runWithTools(AgentRunnerRequest request) {
        // Phase 2: tools XOR structured-output — do NOT force a responseFormat
        // when tools are in play (pass a null schema to the model factory).
        final model = modelFactory.createModel(request.model, DEFAULT_TIMEOUT_SECONDS, null)

        final List<ToolSpecification> specs = request.toolSpecs.collect { ModuleToolAdapter.toToolSpecification(it) }
        final List<ChatMessage> messages = composeMessages(request)
        final int maxIterations = request.maxIterations > 0 ? request.maxIterations : DEFAULT_MAX_ITERATIONS

        log.debug "Running agent model=${request.model}; tools=${specs*.name()}; maxIterations=${maxIterations}"

        AiMessage ai = null
        for( int i = 0; i < maxIterations; i++ ) {
            final ChatRequest req = ChatRequest.builder()
                .messages(messages)
                .toolSpecifications(specs)
                .build()
            final ChatResponse response = model.chat(req)
            ai = response.aiMessage()
            if( !ai.hasToolExecutionRequests() ) {
                // final answer
                return ai.text()
            }
            // append the assistant turn (carrying the tool requests) ...
            messages.add(ai)
            // ... then run each requested tool and append its result
            for( ToolExecutionRequest ter : ai.toolExecutionRequests() ) {
                log.debug "Agent tool call name=${ter.name()}; args=${ter.arguments()}"
                final String result = request.dispatch.call(ter.name(), ter.arguments())
                messages.add(ToolExecutionResultMessage.from(ter, result))
            }
        }

        throw new IllegalStateException("Agent exceeded the maximum number of tool-call iterations (${maxIterations})")
    }

    /**
     * Compose the chat messages: an optional system instruction followed by the
     * user message (the rendered prompt, plus the input record serialized as JSON
     * when present).
     */
    private static List<ChatMessage> composeMessages(AgentRunnerRequest request) {
        String userText = request.prompt
        if( request.inputJson )
            userText += "\n\nInput (JSON):\n" + request.inputJson

        final List<ChatMessage> messages = new ArrayList<ChatMessage>()
        if( request.instruction )
            messages.add(SystemMessage.from(request.instruction))
        messages.add(UserMessage.from(userText))
        return messages
    }
}
