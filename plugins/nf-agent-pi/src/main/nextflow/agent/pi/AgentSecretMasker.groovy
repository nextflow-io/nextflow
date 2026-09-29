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
package nextflow.agent.pi

import java.util.regex.Pattern

import groovy.transform.CompileStatic
import nextflow.SysEnv

/**
 * Masks LLM provider credentials embedded in text captured from an agent runner -- error
 * messages, trace frames, stderr tails -- before any of it reaches the Nextflow log.
 *
 * <p>Masking is keyed off the RESOLVED credential first. A key supplied through
 * {@code agent.apiKey} -- especially spelled {@code secrets.LLM_KEY} -- has no environment
 * variable name, so a name-driven sweep cannot find it. The environment sweep is kept as a
 * backstop for a key the user still exports, and reads through {@link SysEnv} rather than
 * {@code System.getenv} so a test can swap the environment instead of mutating the process.
 * The patterns are a last resort for a credential that is neither: text the provider itself
 * echoed back, or an upstream gateway token.
 *
 * <p>Deliberately shared by both agent transports: an independent copy per transport is exactly
 * how the canonical RPC path came to mask less than the legacy stdin harness. Nothing here is
 * Pi-specific -- it belongs with {@link AgentRpcBroker} if that is ever extracted to a shared
 * plugin.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class AgentSecretMasker {

    /**
     * Provider API-key environment variables swept from captured runner output. The list mirrors
     * the provider tier of the core credential ladder -- {@code AgentConfig} reads exactly these,
     * beside the neutral {@code NXF_AGENT_API_KEY} -- so the redaction backstop and the resolution
     * contract cannot disagree about what counts as a credential.
     *
     * <p>This is the backstop for a key the driver holds but did NOT resolve for the invocation:
     * the one it did resolve is redacted by VALUE, whether it reached the task in band on the RPC
     * start frame or through the {@code env} config scope / {@code secret} directive.
     */
    private static final List<String> SECRET_ENV_KEYS = ['NXF_AGENT_API_KEY','OPENAI_API_KEY','ANTHROPIC_API_KEY','GOOGLE_API_KEY','GEMINI_API_KEY','MISTRAL_API_KEY','OPENROUTER_API_KEY','AZURE_OPENAI_API_KEY']

    private static final Pattern BEARER_TOKEN = ~/(?i)(bearer\s+)[^\s"']+/
    private static final Pattern API_KEY = ~/(?i)\b(sk|rk|pk)-[a-z0-9_-]{8,}/

    /**
     * @param value The text to mask, or {@code null}.
     * @param apiKey The credential resolved for this invocation, or {@code null} when none was.
     * @return The text with every known credential replaced by {@code [REDACTED]}.
     */
    static String redact(String value, String apiKey=null) {
        if( value == null )
            return value
        String result = value
        // The resolved credential first: it is the one value this invocation could actually
        // have leaked, and the only one no name-driven lookup can discover.
        if( apiKey )
            result = result.replace(apiKey, '[REDACTED]')
        for( final key : SECRET_ENV_KEYS ) {
            final secret = SysEnv.get(key)
            if( secret )
                result = result.replace(secret, '[REDACTED]')
        }
        // Mask credentials the text may embed even when the exact forwarded secret is not what
        // leaked (e.g. an upstream Bearer header or a differently-prefixed API key).
        result = result.replaceAll(BEARER_TOKEN, '$1[REDACTED]')
        result = result.replaceAll(API_KEY, '[REDACTED]')
        return result
    }
}
