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

import java.time.Duration

import dev.langchain4j.model.chat.ChatModel
import dev.langchain4j.model.chat.request.ResponseFormat
import dev.langchain4j.model.chat.request.ResponseFormatType
import dev.langchain4j.model.chat.request.json.JsonSchema
import dev.langchain4j.model.openai.OpenAiChatModel
import groovy.transform.CompileStatic

/**
 * Builds a langchain4j {@link ChatModel} from a {@code provider/model} identifier.
 * v1 supports the {@code openai} provider; the API key is read from the
 * provider-standard environment variable.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class ChatModelFactory {

    String apiKey = System.getenv('OPENAI_API_KEY')

    static String providerOf(String modelId) {
        final i = modelId?.indexOf('/') ?: -1
        if( i < 0 )
            throw new IllegalArgumentException("Invalid model id `${modelId}` - expected `provider/model`")
        return modelId.substring(0, i)
    }

    static String modelOf(String modelId) {
        final i = modelId?.indexOf('/') ?: -1
        if( i < 0 )
            throw new IllegalArgumentException("Invalid model id `${modelId}` - expected `provider/model`")
        return modelId.substring(i + 1)
    }

    /**
     * Build a chat model for the given {@code provider/model} id. When a
     * structured-output {@code schema} is provided, the OpenAI model is
     * configured with a strict JSON-schema response format so that the model
     * is constrained to return JSON matching the schema.
     *
     * @param modelId        the {@code provider/model} identifier
     * @param timeoutSeconds the request timeout in seconds
     * @param schema         the structured-output JSON schema, or {@code null}
     *                       for free-form text output
     */
    ChatModel createModel(String modelId, int timeoutSeconds, JsonSchema schema) {
        final provider = providerOf(modelId)
        if( provider != 'openai' )
            throw new IllegalArgumentException("Unsupported agent model provider `${provider}` - v1 supports `openai`")
        if( !apiKey )
            throw new IllegalArgumentException("Missing OPENAI_API_KEY environment variable for agent model `${modelId}`")
        final builder = OpenAiChatModel.builder()
            .apiKey(apiKey)
            .modelName(modelOf(modelId))
            .timeout(Duration.ofSeconds(timeoutSeconds))
        if( schema != null ) {
            final responseFormat = ResponseFormat.builder()
                .type(ResponseFormatType.JSON)
                .jsonSchema(schema)
                .build()
            builder
                .responseFormat(responseFormat)
                .strictJsonSchema(true)
        }
        return builder.build()
    }

    static ChatModel create(String modelId, int timeoutSeconds, JsonSchema schema) {
        new ChatModelFactory().createModel(modelId, timeoutSeconds, schema)
    }
}
