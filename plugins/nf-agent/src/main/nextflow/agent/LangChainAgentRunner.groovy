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

import dev.langchain4j.data.message.ChatMessage
import dev.langchain4j.data.message.SystemMessage
import dev.langchain4j.data.message.UserMessage
import dev.langchain4j.model.chat.ChatModel
import dev.langchain4j.model.chat.request.json.JsonSchema
import dev.langchain4j.model.chat.response.ChatResponse
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import org.pf4j.Extension

/**
 * langchain4j-backed {@link AgentRunner}. v1: single-shot chat (no tool calls).
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@Extension
@CompileStatic
class LangChainAgentRunner implements AgentRunner {

    private static final int DEFAULT_TIMEOUT_SECONDS = 120

    ChatModelFactory modelFactory = new ChatModelFactory()

    @Override
    String run(AgentRunnerRequest request) {
        if( !request.model )
            throw new IllegalArgumentException("Agent `model` directive is required")

        // derive the structured-output schema (when an output record type is declared)
        final JsonSchema schema = request.outputSchema
            ? JsonSchemaMapper.toJsonSchema('Output', request.outputSchema)
            : null

        final model = modelFactory.createModel(request.model, DEFAULT_TIMEOUT_SECONDS, schema)

        // compose the user message: the rendered prompt, followed by the input
        // record serialized as JSON (when present)
        String userText = request.prompt
        if( request.inputJson )
            userText += "\n\nInput (JSON):\n" + request.inputJson

        final List<ChatMessage> messages = new ArrayList<ChatMessage>()
        if( request.instruction )
            messages.add(SystemMessage.from(request.instruction))
        messages.add(UserMessage.from(userText))

        log.debug "Running agent model=${request.model}; messages=${messages.size()}; structured=${schema != null}"
        final ChatResponse response = model.chat(messages)
        return response.aiMessage().text()
    }
}
