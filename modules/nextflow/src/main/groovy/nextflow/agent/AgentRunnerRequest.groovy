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

import groovy.transform.Canonical
import groovy.transform.CompileStatic
import groovy.transform.ToString
import nextflow.exception.AbortOperationException
import nextflow.agent.rpc.AgentRpcHost
import nextflow.agent.rpc.AgentRpcConfig

/**
 * Immutable request passed to an {@link AgentRunner}: the resolved model, the
 * system instruction, the rendered user prompt, the iteration cap, the
 * (currently unused) tool list for forward compatibility, the JSON schema
 * describing the expected structured output, the input record serialized as
 * JSON, the descriptors of the tools the LLM may call, the callback used to
 * execute them, and the optional high-level goal.
 *
 * The {@code toolSpecs} and {@code dispatch} fields are in-JVM only (they carry
 * a live {@link ToolDispatcher} callback) and are never serialized.
 *
 * The agent's tools are PARTITIONED across two fields by who executes them:
 * {@code toolSpecs} carries the BROKERED ones (a descriptor plus a
 * {@code dispatch} callback into the driver), {@code nativeToolNames} the
 * RUNNER-NATIVE ones (bare names the runner serves itself). The two sets are
 * disjoint, and {@link #checkToolPartition} is what makes that structural.
 *
 * Being {@code @Canonical}, the positional constructor order follows the field
 * declaration order below:
 * {@code (model, instruction, prompt, maxIterations, tools, outputSchema, inputJson, toolSpecs, dispatch, requestTimeoutSeconds, goal, agentName, trace, skills, temperature, workDir, apiKey, baseUrl, apiProvider, credentialWithheld, nativeToolNames, brokerHost)}.
 * In practice the request is built with named arguments, so the order is not
 * relied upon at the call site. A NEW field must therefore be APPENDED last, so a
 * positional caller keeps binding the same arguments.
 *
 * The {@code requestTimeoutSeconds} carries the configured per-request LLM chat
 * timeout (from the {@code agent.requestTimeout} config option); when {@code 0}
 * the runner applies its own built-in default.
 *
 * The {@code apiKey} is a CREDENTIAL: it is excluded from the generated
 * {@code toString()} so an interpolated request can never leak it into the log,
 * and it must never be copied into a serialized payload
 * ({@link AgentProtocolSpec#fromRequest} deliberately omits it).
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Canonical
@ToString(excludes='apiKey')
@CompileStatic
class AgentRunnerRequest {
    String model
    String instruction
    String prompt
    int maxIterations
    List tools
    Map outputSchema
    String inputJson
    List<ToolDescriptor> toolSpecs
    ToolDispatcher dispatch
    int requestTimeoutSeconds
    String goal
    /** The agent name, used only to label the execution trace. */
    String agentName
    /** When {@code true}, the runner logs a readable execution trace (turns, model
     * reasoning, tool invocations with inputs/outputs) at INFO level. */
    boolean trace
    /** Portable descriptors of the agent's declared skills. A runner may expose
     * these through its native skill/tool mechanism. Empty/null when no
     * {@code skills} directive is declared. */
    List<SkillDescriptor> skills
    /** Sampling temperature applied to the LLM chat, or {@code null} to leave the
     * provider default. Left UNSET (null) on BOTH the tool-free task path and the
     * legacy tool/skill path; an explicit value is applied only when a caller opts in
     * ({@link ChatModelFactory} calls {@code .temperature(...)} only when non-null).
     * It is deliberately NOT pinned to {@code 0.0}: reasoning models (e.g. gpt-5-mini,
     * used by every committed example) reject an explicit temperature with HTTP 400,
     * and resume is pure input-keyed memoization (a hit replays the stored generation
     * without calling the model), so replay correctness does not depend on temperature=0.
     * Declared LAST to preserve the positional {@code @Canonical} constructor order
     * used by existing callers (a trailing omitted arg defaults to {@code null}). */
    Double temperature
    /** Portable string form of the invocation task work directory. */
    String workDir
    /** The LLM provider credential resolved by the core ladder and scoped to this
     * model's provider ({@link AgentConfig#apiKeyFor}), or {@code null} when nothing
     * resolves. A runner MUST NOT read the environment itself, and MUST present
     * {@link #credential()} rather than this field. Declared LAST-BUT-ONE to
     * preserve the positional {@code @Canonical} constructor order. */
    String apiKey
    /** The OpenAI-compatible endpoint resolved by the core ladder
     * ({@link AgentConfig#baseUrlFor}), or {@code null} to leave the provider
     * default. Not a secret: unlike {@code apiKey} it travels on the protocol spec
     * and enters the resume cache key. */
    String baseUrl
    /** The provider NAMESPACE the core ladder resolved the pair above from
     * ({@link AgentConfig#apiProviderFor}) -- {@code agent.apiProvider}, else the provider of a
     * well-known endpoint host, else the model-id prefix. Not the wire protocol, and not a secret.
     * Carried so a runner that must report a missing credential can name the variables the ladder
     * ACTUALLY consulted instead of asserting the OpenAI ones. */
    String apiProvider
    /** {@code true} when a {@code <PROVIDER>_API_KEY} DID resolve for this model and the endpoint
     * gate withheld it ({@link AgentConfig#credentialWithheldFor}); {@code false} both when a
     * credential travels and when none exists at all.
     *
     * <p>It is what separates a diagnosable misconfiguration from the ordinary no-credential case,
     * and {@link #credential()} refuses to paper over it with {@link #PLACEHOLDER_API_KEY}. */
    boolean credentialWithheld
    /** Wire names of the RUNNER-NATIVE tools the agent selected -- the {@code fs:} and
     * {@code shell:} leaves of the tool grammar, as the bare names the model sees
     * ({@code read}, {@code write}, {@code edit}, {@code ls}, {@code grep}, {@code find},
     * {@code bash}).
     *
     * <p>They travel BESIDE {@code toolSpecs} and never inside it. A native tool is served by the
     * RUNNER -- a pi SDK builtin enabled through the session allowlist, an in-JVM implementation on
     * langchain4j -- so it has no descriptor, no schema of ours and no dispatcher, and there is
     * nothing for the driver to execute on its behalf. The separation is the safety property: the
     * broker authorizes exactly the brokered names ({@link #brokeredToolNames}), so a name that
     * stays out of {@code toolSpecs} cannot be called back into the driver JVM at all.
     *
     * <p>A runner that serves none of them (there is no {@code shell:} on langchain4j) simply
     * ignores the field; core has already refused the refs that runner cannot honour.
     *
     * <p>Declared LAST to preserve the positional {@code @Canonical} constructor order. */
    List<String> nativeToolNames

    /** The address a CONTAINERIZED agent task must dial to reach the driver's RPC broker, carried
     * with the ladder row that produced it ({@link AgentRpcHost}).
     *
     * <p>It rides on the request because it is resolved PER AGENT DEFINITION, pre-ignition, from
     * context the broker does not have -- the executor instance, the engine config, the task's
     * container options. A run may legitimately hold several: one agent on the local docker engine
     * resolves the engine host alias while another on {@code k8s} resolves the driver's pod address,
     * and each task must be told ITS OWN or it dials an address it cannot route to and the driver
     * holds its capability for the full {@code agent.rpc.capabilityTimeout}.
     *
     * <p>{@code null} on the in-JVM runner path, where nothing dials anything, and on a runner that
     * registers without the guard -- the broker falls back to
     * {@link AgentRpcConfig#resolveBrokerHost} there. Never part of the resume cache key: it is a
     * property of the machine, not of the agent. Declared LAST to preserve the positional
     * {@code @Canonical} constructor order. */
    AgentRpcHost brokerHost

    /**
     * Stand-in credential presented when an endpoint is declared but no credential resolved
     * (design D8). A local vLLM or Ollama needs none, yet both consumers require SOMETHING: the
     * langchain4j OpenAI client rejects an empty key, and the pi runner fails the run outright
     * (verified: {@code No API key found for openai} from the session, or
     * {@code Provider is not configured} from {@code prepareRequest}) -- which is exactly the
     * local-first case this feature exists to unblock. Cosmetic on langchain4j, load-bearing on pi.
     */
    static final String PLACEHOLDER_API_KEY = 'nxf-no-credential'

    /**
     * Fail unless the brokered and the runner-native halves are DISJOINT.
     *
     * <p>This is enforcement, not a sanity check. Both halves come from one resolved selection, so
     * an overlap is a bug in whoever assembled the request -- but it is the one bug whose
     * consequence is a boundary crossing rather than a wrong answer: a native name reaching
     * {@code toolSpecs} enters the broker's authorization set, and a tool the model believes runs
     * inside the runner container would then be executed in the DRIVER JVM. Better a run that
     * aborts before a job is submitted than one that silently relocates {@code bash}.
     *
     * <p>Called from the two chokepoints every runner passes through -- the portable payload
     * builder ({@link AgentProtocolSpec#fromRequest}) and the broker registration that mints the
     * allowlist -- so the invariant holds by construction rather than by convention.
     */
    void checkToolPartition() {
        if( !nativeToolNames || !toolSpecs )
            return
        final Set<String> brokered = brokeredNames()
        final List<String> overlap = nativeToolNames.findAll { String it -> brokered.contains(it) }
        if( overlap )
            throw new IllegalStateException("Agent tool partition violated: ${overlap.join(', ')} - a runner-native tool is served by the runner itself and must never be carried as a brokered tool descriptor")
    }

    /**
     * The tool names a runner is authorized to call BACK into the driver with -- the brokered ones,
     * and only those. Validated: see {@link #checkToolPartition}.
     */
    Set<String> brokeredToolNames() {
        checkToolPartition()
        return brokeredNames()
    }

    private Set<String> brokeredNames() {
        final Set<String> names = new LinkedHashSet<String>()
        for( ToolDescriptor it : (toolSpecs ?: Collections.<ToolDescriptor>emptyList()) )
            names.add(it.name)
        return names
    }

    /**
     * The credential the runner must actually present for this request: the resolved one, or the
     * D8 placeholder when an endpoint is declared and nothing resolved, or {@code null} when
     * neither is set.
     *
     * <p>The placeholder is confined to the OpenAI protocol on purpose. A runner installs what it
     * is given as the credential OF THE MODEL'S PROVIDER and that ownership beats the ambient
     * environment (see {@link AgentConfig#apiKeyFor}), so a placeholder sent for, say,
     * {@code anthropic/claude-*} would MASK an exported {@code ANTHROPIC_API_KEY} the runner can
     * resolve by itself, turning a working run into a 401. An openai-protocol endpoint is also the
     * only shape the no-credential case has: it is a local vLLM/Ollama/llama.cpp server.
     *
     * <p>{@code null} is NOT an error here: a runner may have credential sources core cannot see
     * (pi reads its own store and provider-specific variables), so rejecting the request is the
     * runner's decision, not this class'.
     */
    String credential() {
        // A WITHHELD credential is not a missing one. The endpoint gate refused to send a key that
        // resolved (AgentConfig.credentialWithheldFor), which is a misconfiguration a runner can
        // name precisely; substituting the placeholder here would convert that diagnosis into the
        // opaque 401 the placeholder rule exists to avoid, and on the pi path it would additionally
        // shadow whatever the container was given out of band. The placeholder is ONLY for a
        // genuine no-credential local endpoint.
        if( credentialWithheld )
            return null
        return AgentConfig.isOpenAiProtocol(model)
            ? credentialFor(apiKey, baseUrl)
            : apiKey
    }

    /**
     * The same rule for a caller that holds the resolved pair but not a request -- the langchain4j
     * {@code ChatModelFactory}, which has already established that the model is OpenAI-protocol.
     * Kept here so the two runners cannot drift into two placeholder rules.
     *
     * <p>This is the CORE-IS-THE-ONLY-SOURCE path: it answers "what must I present", so it is also
     * where design D5 lives. The placeholder assumes the endpoint needs no credential, which is
     * true of a local vLLM/Ollama and false of a provider's own API, so a {@code baseUrl} whose
     * host is a KNOWN provider endpoint fails here with a diagnosis instead of buying an opaque
     * 401 one request later. A runner with credential sources of its own must read
     * {@link #apiKey} and never reach this: core resolving nothing is not an error for it (a
     * containerized pi task may be given the key by {@code env}/{@code secret} or a Kubernetes
     * Secret the driver cannot see), so this must not become a run-aborting check at agent build
     * time.
     */
    static String credentialFor(String apiKey, String baseUrl) {
        if( apiKey )
            return apiKey
        if( !baseUrl )
            return null
        final provider = AgentConfig.inferProviderFromUrl(baseUrl)
        if( provider )
            // D5 names all three ways out. `agent.apiProvider` is the third because the namespace
            // this message reports was INFERRED from the endpoint host, and an inference is exactly
            // the kind of answer a user may need to override.
            throw new AbortOperationException("Missing `${provider}` credential for the agent endpoint ${baseUrl} - ${AgentConfig.missingCredentialHint(provider)}; the `${provider}` namespace was inferred from the endpoint host, so set `agent.apiProvider` to name a different one")
        return PLACEHOLDER_API_KEY
    }
}
