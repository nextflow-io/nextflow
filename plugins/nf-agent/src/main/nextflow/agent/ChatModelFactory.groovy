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
import dev.langchain4j.model.chat.listener.ChatModelListener
import dev.langchain4j.model.chat.request.ResponseFormat
import dev.langchain4j.model.chat.request.ResponseFormatType
import dev.langchain4j.model.chat.request.json.JsonSchema
import dev.langchain4j.model.openai.OpenAiChatModel
import groovy.transform.CompileStatic

/**
 * Builds a langchain4j {@link ChatModel} from a {@code provider/model} identifier.
 *
 * The {@code openai} prefix names the OpenAI WIRE PROTOCOL, not the vendor: with an
 * {@code agent.baseUrl} it also reaches any OpenAI-compatible endpoint (vLLM, Ollama,
 * llama.cpp, OpenRouter, LiteLLM, a corporate gateway).
 *
 * This factory NEVER reads the environment. The endpoint and the credential are
 * resolved once in core ({@code AgentConfig.apiKeyFor}/{@code baseUrlFor} implement the
 * config -> {@code NXF_AGENT_*} -> {@code <PROVIDER>_*} ladder, where the provider namespace
 * is {@code agent.apiProvider}, else the provider of a well-known endpoint host, else the
 * model-id prefix), travel on the {@link nextflow.agent.AgentRunnerRequest} and are passed in
 * as {@code createModel} arguments — so one ladder is authoritative for every runner. The
 * no-credential rule ({@code AgentRunnerRequest.PLACEHOLDER_API_KEY}) is shared with the pi
 * runner for the same reason, and with it the two rules that bound it: a {@code baseUrl} whose
 * host IS a known provider's endpoint fails here rather than presenting a placeholder that
 * buys an opaque 401, and a credential that RESOLVED and was withheld by the endpoint gate
 * ({@code AgentConfig.credentialWithheldFor}) fails too — the placeholder is for a genuine
 * no-credential local endpoint, never for a misroute. This runner has no credential source of
 * its own, so both are fatal here; the pi runner only warns.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class ChatModelFactory {

    private static int slashIndex(String modelId) {
        final i = modelId?.indexOf('/') ?: -1
        if( i < 0 )
            throw new IllegalArgumentException("Invalid model id `${modelId}` - expected `provider/model`")
        return i
    }

    /**
     * The provider prefix, LOWER-CASED so this and {@link AgentConfig#providerPrefixOf} give the
     * same answer for the same model id. They must: core resolves and TRANSMITS a credential for
     * the provider it computes, so a mixed-case {@code OpenAI/gpt-4o} that core reads as
     * {@code openai} and this factory read as {@code OpenAI} disagreed about whose key was in
     * flight -- core sent one, and the runner rejected the model as an unsupported provider.
     */
    static String providerOf(String modelId) {
        return modelId.substring(0, slashIndex(modelId)).toLowerCase()
    }

    static String modelOf(String modelId) {
        return modelId.substring(slashIndex(modelId) + 1)
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
     * @param temperature    the sampling temperature to apply, or {@code null} to
     *                       leave the provider default (legacy tool/skill path)
     * @param apiKey         the credential resolved by core, or {@code null} when none
     *                       resolved (allowed only together with a {@code baseUrl})
     * @param baseUrl        the OpenAI-compatible endpoint resolved by core, or
     *                       {@code null} to use the provider default
     * @param apiProvider    the provider NAMESPACE core resolved the pair from
     *                       ({@code AgentRunnerRequest.apiProvider}), used to name the right
     *                       variables in a missing-credential message; {@code null} when unknown
     * @param credentialWithheld  {@code true} when a provider credential resolved and the endpoint
     *                       gate refused to send it ({@code AgentConfig.credentialWithheldFor})
     * @param listeners      optional chat-model listeners (e.g. the execution
     *                       tracer), or {@code null} for none
     *
     * NOTE: the credential parameters are declared AFTER {@code temperature} and BEFORE the
     * defaulted {@code listeners} parameter, so the generated default-argument overload cannot
     * bind one of them into the listener slot.
     */
    ChatModel createModel(String modelId, int timeoutSeconds, JsonSchema schema, Double temperature, String apiKey, String baseUrl, String apiProvider, boolean credentialWithheld, List<ChatModelListener> listeners = null) {
        final provider = providerOf(modelId)
        if( provider != AgentConfig.OPENAI_PROVIDER )
            throw new IllegalArgumentException("Unsupported agent model provider `${provider}` - only the `openai` wire protocol is supported; to reach an OpenAI-compatible endpoint keep the `openai/` prefix and set `agent.baseUrl`")
        // The credential namespace the ladder consulted, as core resolved it. Falling back to the
        // prefix keeps a request built before this field existed (or by a test) from claiming a
        // namespace nobody chose.
        final namespace = apiProvider ?: provider
        if( credentialWithheld )
            // A key RESOLVED and the endpoint gate refused to send it. That is a misconfiguration
            // with a name, so say it here rather than let the placeholder below turn it into a 401.
            // langchain4j has no other credential source, so this is fatal -- deliberately UNLIKE
            // the pi path, which only warns (see AgentRpcBroker.warnWithheldCredential).
            throw new IllegalArgumentException("The `${namespace}` credential resolved for agent model `${modelId}` was withheld from ${baseUrl ? "the endpoint ${baseUrl}" : "the default `${provider}` endpoint"} - that endpoint is not one `${namespace}` owns; set `agent.apiKey` in the Nextflow configuration, or the `NXF_AGENT_API_KEY` environment variable, to the credential it does accept, or `agent.apiProvider` to name the namespace it belongs to")
        // the shared D8 rule: the resolved credential, else a placeholder when an endpoint is
        // declared (the endpoint is assumed to need none: a local vLLM/Ollama), else nothing --
        // and, since the endpoint may instead be a known provider's own API, an abort when that
        // assumption is plainly false. Owned by AgentRunnerRequest so this runner and the pi
        // runner cannot drift apart.
        final credential = AgentRunnerRequest.credentialFor(apiKey, baseUrl)
        if( !credential )
            // Name the variables the ladder ACTUALLY consulted for the resolved namespace. Asserting
            // OPENAI_API_KEY here was wrong whenever `agent.apiProvider` redirected the namespace --
            // the very thing the old comment claimed to avoid. The namespace is still named as a
            // remedy, because it is also the rung a user may need to CHANGE.
            throw new IllegalArgumentException("Missing LLM provider credential for agent model `${modelId}` - ${AgentConfig.missingCredentialHint(namespace)}; the credential namespace is `${namespace}`, set `agent.apiProvider` to name a different one")
        final builder = OpenAiChatModel.builder()
            .apiKey(credential)
            .modelName(modelOf(modelId))
            .timeout(Duration.ofSeconds(timeoutSeconds))
        if( baseUrl )
            builder.baseUrl(baseUrl)
        if( temperature != null )
            builder.temperature(temperature)
        if( listeners ) {
            builder.listeners(listeners)
            // also request the model's reasoning content: it is parsed into AiMessage.thinking()
            // for models/providers that return it (reasoning models, or providers exposing
            // `reasoning_content`), and is absent otherwise. Only returnThinking is set —
            // sendThinking stays false, so the reasoning is never echoed back to the model on
            // later turns (no risk of chat-completions rejecting an echoed reasoning field).
            builder.returnThinking(true)
        }
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
}
