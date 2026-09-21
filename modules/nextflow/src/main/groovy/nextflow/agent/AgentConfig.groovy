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

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.SysEnv
import nextflow.config.spec.ConfigOption
import nextflow.config.spec.ConfigScope
import nextflow.config.spec.ScopeName
import nextflow.config.spec.SpecNode
import nextflow.exception.AbortOperationException
import nextflow.script.dsl.ConfigSelectorResolver
import nextflow.script.dsl.Description
import nextflow.script.dsl.ProcessConfigBuilder
import nextflow.util.Duration
import nextflow.util.MemoryUnit
import nextflow.agent.rpc.AgentRpcConfig

/**
 * Model the agent-only options of the `agent` configuration scope: the runner, the endpoint and
 * credential of the LLM provider, and the model/iteration/timeout/trace settings that are
 * DEFAULTS for agent declarations.
 *
 * <p>The rest of the scope is the task-directive axis (the same directives a process
 * accepts, e.g. {@code executor}, {@code cpus}, {@code container}) applied to the agent
 * task by {@link nextflow.script.dsl.ProcessConfigBuilder}, independently from the
 * {@code process} scope. Both axes are resolved with the same selector ladder — see
 * {@link #resolveOptions} and {@link nextflow.script.AgentDef#buildAgentTask}.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@ScopeName("agent")
@Description("""
    The `agent` scope controls the default settings applied to agents.
""")
@Slf4j
@CompileStatic
class AgentConfig implements ConfigScope {

    static final String DEFAULT_EXECUTOR = 'local'

    /**
     * Default cap on concurrently running agent tasks, applied when the agent declares no
     * {@code maxForks} of its own.
     *
     * <p>Unlike a process, an agent is admitted by {@link nextflow.processor.AgentPollingMonitor},
     * which applies no cpu or capacity throttle -- so without this an agent fans out as wide as its
     * input channel. That is a THROUGHPUT concern rather than a correctness one (orchestration runs
     * on its own pool, see {@link nextflow.Session#getAgentExecService}), but an uncapped fan-out
     * of concurrent LLM calls invites provider rate-limiting.
     *
     * <p>Deliberately a constant rather than a function of the core count: an agent waits on a
     * remote model, so the machine's cpus say nothing about how many calls a provider tolerates.
     * Conservative by intent and meant to be tuned -- raise it with {@code maxForks} on the agent,
     * or in the {@code agent} scope.
     */
    static final int DEFAULT_MAX_FORKS = 10

    /**
     * The provider prefix that names the OpenAI WIRE PROTOCOL (design D1). It no longer gates the
     * endpoint/credential ladder -- {@link #apiProviderFor} answers "whose key and whose endpoint"
     * for every provider -- but it still answers "which client speaks to it", the one question the
     * model-id prefix is authoritative for.
     */
    static final String OPENAI_PROVIDER = 'openai'

    /**
     * The provider tokens the PROVIDER tier of the ladder knows, each mapped to the environment
     * variables read for it, in order -- first TRUTHY hit wins (design D2). A provider is not
     * always 1:1 with a variable ({@code GOOGLE_API_KEY} / {@code GEMINI_API_KEY}), hence a list.
     *
     * <p>The namespace is CLOSED on purpose. The provider is pipeline-supplied text (the
     * {@code agent.apiProvider} option, or the model-id prefix) and the endpoint that would receive
     * the credential is pipeline-supplied config too, so uppercasing an arbitrary token into
     * {@code <X>_API_KEY} would let a pipeline name ANY variable in the driver's environment and
     * have it delivered to an address of its choosing. An unlisted token resolves nothing; written
     * as {@code agent.apiProvider} it is rejected outright (see {@link #resolveApiProvider}).
     *
     * <p>The variables match {@code AgentSecretMasker.SECRET_ENV_KEYS}, which already redacts
     * exactly these from captured runner output: the redaction backstop and the resolution contract
     * must not disagree about what counts as a credential.
     */
    private static final Map<String,List<String>> PROVIDER_API_KEY_VARS = Collections.unmodifiableMap([
            anthropic : List.of('ANTHROPIC_API_KEY'),
            azure     : List.of('AZURE_OPENAI_API_KEY'),
            gemini    : List.of('GEMINI_API_KEY', 'GOOGLE_API_KEY'),
            google    : List.of('GOOGLE_API_KEY', 'GEMINI_API_KEY'),
            mistral   : List.of('MISTRAL_API_KEY'),
            openai    : List.of('OPENAI_API_KEY'),
            openrouter: List.of('OPENROUTER_API_KEY') ] as Map<String,List<String>>)

    /**
     * The endpoint half of {@link #PROVIDER_API_KEY_VARS}, deliberately NOT the same key set:
     * {@code OPENAI_BASE_URL} and {@code ANTHROPIC_BASE_URL} are the official SDK spellings and
     * {@code AZURE_OPENAI_ENDPOINT} is Azure's, while Google, Mistral and OpenRouter define no such
     * variable at all. Reading a {@code MISTRAL_BASE_URL} would be Nextflow inventing a convention
     * under a vendor's name and then owing it forever; retargeting one of those providers is what
     * {@code agent.baseUrl} and {@code NXF_AGENT_BASE_URL} are for.
     */
    private static final Map<String,List<String>> PROVIDER_BASE_URL_VARS = Collections.unmodifiableMap([
            anthropic: List.of('ANTHROPIC_BASE_URL'),
            azure    : List.of('AZURE_OPENAI_ENDPOINT'),
            openai   : List.of('OPENAI_BASE_URL') ] as Map<String,List<String>>)

    /**
     * Hosts whose provider is unambiguous, used to infer the credential namespace from a
     * provider-neutral endpoint (design D3). Deliberately short: every entry is a compatibility
     * commitment, and a wrong one misroutes a credential. Matched on the HOST only -- exactly or as
     * a dot-suffix -- never as a substring, so neither {@code https://evil.example/openai/v1} nor
     * {@code https://api.openai.com.evil.example} matches.
     */
    private static final Map<String,String> PROVIDER_HOSTS = Collections.unmodifiableMap([
            'api.openai.com'   : 'openai',
            'api.anthropic.com': 'anthropic',
            'openrouter.ai'    : 'openrouter',
            'api.mistral.ai'   : 'mistral' ] as Map<String,String>)

    /**
     * The `agent` scope options that are NOT task directives, derived from the
     * {@code @ConfigOption} fields and the nested {@link ConfigScope} fields declared by this
     * class. Task directives come from {@link nextflow.script.dsl.ProcessDsl.DirectiveDsl}
     * instead, so a new agent-only option can never be mistaken for a process directive.
     *
     * NOTE: a nested {@code agent.<x>} scope is declared as a field whose type implements
     * {@link ConfigScope} and which is NOT annotated {@code @ConfigOption} (see {@code rpc}) --
     * {@link nextflow.config.spec.SpecNode.Scope#of} recurses into it, so every option inside is
     * individually validated. Such a field MUST also be picked up here, otherwise
     * {@code agent { x { ... } }} is applied as a bogus task directive.
     *
     * A plugin CANNOT contribute a nested scope: {@link nextflow.config.ConfigValidator} registers a
     * {@code @ScopeName} as a single top-level key, so a dotted {@code @ScopeName('agent.x')} is
     * unreachable for the segment-by-segment spec lookup. Nested agent scopes therefore live here.
     */
    static final Set<String> AGENT_ONLY_OPTIONS = declaredOptions()

    private static Set<String> declaredOptions() {
        final names = SpecNode.Scope.of(AgentConfig, '').children().keySet()
        return Collections.unmodifiableSet(new LinkedHashSet<String>(names))
    }

    @ConfigOption
    @Description("""
        The agent runner extension to use (for example `pi` or `langchain4j`). Required when more than one runner plugin is enabled.
    """)
    final String runner

    @ConfigOption
    @Description("""
        The default model (`provider/model`) used by an agent that does not declare a `model` directive.
    """)
    final String model

    @ConfigOption
    @Description("""
        The provider namespace the API key and base URL are taken from when they come from the environment, one of `anthropic`, `azure`, `gemini`, `google`, `mistral`, `openai`, `openrouter`. When not specified, it is inferred from the host of `baseUrl` when that names a well-known provider, and otherwise from the model id prefix. It does NOT select the wire protocol, which is always the model id prefix.
    """)
    final String apiProvider

    @ConfigOption
    @Description("""
        The API key used to authenticate with the LLM provider. When not specified, the `NXF_AGENT_API_KEY` environment variable is used, then the API provider's own variable (e.g. `OPENAI_API_KEY`, `ANTHROPIC_API_KEY`) when the endpoint it would be sent to belongs to that provider.
    """)
    final String apiKey

    @ConfigOption
    @Description("""
        The base URL of the endpoint serving the model, e.g. a local vLLM or Ollama server or a gateway. When not specified, the `NXF_AGENT_BASE_URL` environment variable is used, then the API provider's own variable (`OPENAI_BASE_URL`, `ANTHROPIC_BASE_URL`, `AZURE_OPENAI_ENDPOINT`), then the provider default.
    """)
    final String baseUrl

    @ConfigOption
    @Description("""
        The default maximum number of LLM iterations used by an agent that does not declare a `maxIterations` directive (default: `20`).
    """)
    final Integer maxIterations

    @ConfigOption
    @Description("""
        The amount of time to wait for an LLM chat request to complete before failing (default: `120 sec`).
    """)
    final Duration requestTimeout

    @ConfigOption
    @Description("""
        The maximum size of a structured tool-output file whose contents are passed to the LLM; larger outputs are returned as a path handle (default: `32 KB`).
    """)
    final MemoryUnit maxToolOutputInlineSize

    @ConfigOption
    @Description("""
        When `true`, log a readable trace of each agent's execution - the turns, the model reasoning and the tool invocations - at INFO level; the tool inputs and outputs are logged at DEBUG. Enabled by the `-with-agent-trace` run option (default: `false`).
    """)
    final Boolean trace

    @Description("""
        Settings for the agent RPC broker that a canonical agent task connects back to.
    """)
    final AgentRpcConfig rpc

    /**
     * The environment the PROVIDER tier of the ladder reads. Held rather than re-read because that
     * tier is resolved per MODEL, long after construction, and an object whose answers depend on
     * when they are asked cannot be reasoned about (or swapped in a test).
     *
     * NOTE: must stay un-annotated and NOT a {@link ConfigScope}, or {@link #declaredOptions}
     * would publish it as a bogus {@code agent} option.
     */
    private final Map sysEnv

    /**
     * The provider inferred from the provider-NEUTRAL endpoint (design D3), or {@code null} when
     * there is none or its host is not recognized. Eager because it depends only on
     * {@link #baseUrl}; the third rung of the provider ladder (the model prefix) cannot be, see
     * {@link #apiProviderFor}.
     */
    private final String inferredProvider

    /* required by extension point -- do not remove */
    AgentConfig() {}

    AgentConfig(Map opts) {
        this(opts, SysEnv.get())
    }

    /**
     * The same, with the environment passed in rather than read: a test swaps a map instead of the
     * process environment, and the deferred provider tier is guaranteed to see the SAME environment
     * the neutral tiers were resolved from.
     */
    AgentConfig(Map opts, Map env) {
        sysEnv = env
        runner = opts.runner as String
        model = opts.model as String
        // The endpoint and the credential resolve on ONE ladder for every provider (design D2):
        // the config option, then the provider-neutral `NXF_AGENT_*` variable, then the provider's
        // own `<PROVIDER>_*` variable. Only the first two can be resolved HERE -- the third needs
        // the provider, which needs the effective model, which this class never sees (the `model`
        // directive is joined with the `agent.model` default one layer up, in AgentDef). So these
        // fields carry the NEUTRAL tiers and the provider tier lives in #apiKeyFor/#baseUrlFor,
        // the accessors a runner is actually handed.
        apiKey = resolveNeutralApiKey(env, opts)
        baseUrl = resolveNeutralBaseUrl(env, opts)
        apiProvider = resolveApiProvider(opts)
        // inference reads the NEUTRAL endpoint only: `<PROVIDER>_BASE_URL` cannot feed the
        // inference that decides which provider's variable to read in the first place (D3)
        inferredProvider = inferProviderFromUrl(baseUrl)
        // D1 lets an explicit provider win, but one pointed at ANOTHER provider's own host is far
        // likelier a mistake than an intention -- and the mistake ships a credential to a third
        // party. Say so here, once per agent; #apiKeyFor withholds the provider tier in that case.
        if( apiProvider && inferredProvider && apiProvider != inferredProvider )
            log.warn "Agent config sets `agent.apiProvider = '${apiProvider}'` but the endpoint ${baseUrl} is a known `${inferredProvider}` endpoint - no `${apiProvider}` credential will be sent to it; set `agent.apiKey` to the credential this endpoint accepts"
        // null casts to null for object types, so the coercion doubles as the null-guard
        maxIterations = opts.maxIterations as Integer
        requestTimeout = opts.requestTimeout as Duration
        maxToolOutputInlineSize = opts.maxToolOutputInlineSize as MemoryUnit
        trace = opts.trace as Boolean
        rpc = new AgentRpcConfig(opts.rpc instanceof Map ? (Map)opts.rpc : Collections.emptyMap())
    }

    /**
     * Resolve the agent-only options for one agent, applying the same selector ladder as
     * {@link nextflow.script.dsl.ProcessConfigBuilder#applyConfig}: plain scope (weakest), then
     * `withLabel:`, then `withName:` on the base name, the alias, and the fully-qualified name
     * (strongest). Last write wins; unlike the directive axis there is no `ext` merge, no
     * repeatable directive and no built-in default to preserve, so a flat overwrite is enough.
     *
     * @param scope         the raw `agent` config scope
     * @param labels        the agent's declared labels (may be null/empty)
     * @param baseName      the agent's declared name
     * @param simpleName    the agent's (possibly aliased) simple name
     * @param fqName        the agent's fully-qualified name
     */
    static Map<String,Object> resolveOptions(Map<String,?> scope, List<String> labels,
            String baseName, String simpleName, String fqName) {
        final result = new LinkedHashMap<String,Object>()
        copyOptions(scope, result)
        final targets = labels ?: List.of('')
        for( final settings : ConfigSelectorResolver.matchingSettings(scope, targets, baseName, simpleName, fqName) )
            copyOptions(settings, result)
        return result
    }

    /**
     * Copy the agent-only options out of a plain scope or a selector body. A non-map selector
     * body is skipped, NOT cast: this method runs before {@code applyConfig}, which reports it
     * with the canonical {@code ConfigParseException} ("Unknown config settings for agent with
     * name: ...").
     */
    private static void copyOptions(Object source, Map<String,Object> target) {
        if( !(source instanceof Map) )
            return
        for( final entry : ((Map<?,?>) source).entrySet() ) {
            final key = entry.key.toString()
            if( AGENT_ONLY_OPTIONS.contains(key) )
                target.put(key, entry.value)
        }
    }

    /**
     * The credential tiers that are PROVIDER-NEUTRAL: the {@code agent.apiKey} option
     * (selector-aware, already merged into {@code opts} by {@link #resolveOptions}) and
     * {@code NXF_AGENT_API_KEY}. Both are named by the user for THIS agent's provider, whichever
     * that is, so they may be handed to any runner and to any endpoint; the provider tier may not.
     * See {@link #apiKeyFor}.
     *
     * <p>Returns {@code null} when nothing resolves: whether a credential is REQUIRED is a runner
     * concern (a local endpoint needs none), so this must never throw -- it is reached from the
     * constructor for every agent, including those driven by a test runner.
     *
     * <p>Every tier is tested for TRUTHINESS, not for null, so an empty value falls through
     * instead of shadowing the tiers below it: {@code export NXF_AGENT_API_KEY=} and
     * {@code agent.apiKey = params.key} with {@code params.key} unset both yield {@code ''}, and
     * neither is a credential. Same rule as the config tier of {@code AwsConfig.getAwsProfile0}.
     * The user-visible consequence: an explicitly EMPTY {@code agent.apiKey} does not mean "this
     * endpoint needs no credential" -- it falls through to the environment. The no-credential path
     * of design D8 is reached by leaving the option UNSET while setting a {@code baseUrl} that is
     * not a known provider's host (design D5, {@link AgentRunnerRequest#credentialFor}).
     *
     * @param env  the environment to read, always {@link SysEnv#get()} in production
     * @param opts the resolved agent-only options
     */
    static protected String resolveNeutralApiKey(Map env, Map opts) {
        return (opts?.apiKey as String)
            ?: (env?.get('NXF_AGENT_API_KEY') as String)
            ?: null
    }

    /** The endpoint tiers that are PROVIDER-NEUTRAL -- see {@link #resolveNeutralApiKey}. */
    static protected String resolveNeutralBaseUrl(Map env, Map opts) {
        return (opts?.baseUrl as String)
            ?: (env?.get('NXF_AGENT_BASE_URL') as String)
            ?: null
    }

    /**
     * The {@code agent.apiProvider} option, normalized to the token the provider tables are keyed
     * by, or {@code null} when unset (empty is unset, by the same truthiness rule as every other
     * tier).
     *
     * <p>An unknown value is REJECTED rather than uppercased into a variable name: the closed
     * namespace is what keeps {@code <PROVIDER>_API_KEY} from naming arbitrary environment
     * variables (see {@link #PROVIDER_API_KEY_VARS}), and an unrecognized token would otherwise
     * land as a silent tier-3 miss -- the failure mode a typo produces, diagnosed nowhere.
     */
    static protected String resolveApiProvider(Map opts) {
        final value = (opts?.apiProvider as String)?.trim()?.toLowerCase()
        if( !value )
            return null
        if( !isKnownProvider(value) )
            throw new AbortOperationException("Invalid `agent.apiProvider` value `${opts.apiProvider}` - expected one of: ${knownProviders().join(', ')}")
        return value
    }

    /** The provider tokens the {@code <PROVIDER>_*} tier knows. */
    static Set<String> knownProviders() {
        return PROVIDER_API_KEY_VARS.keySet()
    }

    static boolean isKnownProvider(String provider) {
        return provider != null && PROVIDER_API_KEY_VARS.containsKey(provider)
    }

    /** The ordered {@code <PROVIDER>_API_KEY} candidates of a provider; empty when unknown. */
    static List<String> apiKeyVarsFor(String provider) {
        return PROVIDER_API_KEY_VARS.getOrDefault(provider, Collections.<String>emptyList())
    }

    /** The ordered endpoint candidates of a provider; empty when it defines no such variable. */
    static List<String> baseUrlVarsFor(String provider) {
        return PROVIDER_BASE_URL_VARS.getOrDefault(provider, Collections.<String>emptyList())
    }

    /** The name of the first {@code <PROVIDER>_API_KEY} candidate SET in {@code env}, or null. */
    static protected String providerApiKeyVar(Map env, String provider) {
        return firstSetVar(env, apiKeyVarsFor(provider))
    }

    /** The PROVIDER tier of the credential ladder, ungated -- {@link #apiKeyFor} owns the gate. */
    static protected String resolveProviderApiKey(Map env, String provider) {
        final name = providerApiKeyVar(env, provider)
        return name != null ? env.get(name) as String : null
    }

    /** The PROVIDER tier of the endpoint ladder. */
    static protected String resolveProviderBaseUrl(Map env, String provider) {
        final name = firstSetVar(env, baseUrlVarsFor(provider))
        return name != null ? env.get(name) as String : null
    }

    private static String firstSetVar(Map env, List<String> names) {
        if( env == null )
            return null
        for( final name : names )
            if( env.get(name) )     // TRUTHINESS, so an exported-but-empty variable falls through
                return name
        return null
    }

    /**
     * The provider that owns the given endpoint, or {@code null} when its host is not one this
     * table recognizes (design D3). Never guesses and never partial-matches: no match means the
     * caller falls back to the model prefix, which is the conservative answer.
     */
    static String inferProviderFromUrl(String url) {
        final host = hostOf(url)
        if( !host )
            return null
        for( final entry : PROVIDER_HOSTS.entrySet() ) {
            if( host == entry.key || host.endsWith('.' + entry.key) )
                return entry.value
        }
        return null
    }

    /**
     * The lower-case host of a URL, or {@code null} when it has none (a relative or malformed
     * value). A trailing dot -- the absolute-FQDN spelling {@code api.openai.com.} -- is stripped
     * so it cannot slip past an exact match.
     */
    private static String hostOf(String url) {
        if( !url )
            return null
        try {
            String host = new URI(url.trim()).getHost()
            if( host == null )
                return null
            host = host.toLowerCase()
            return host.endsWith('.') ? host.substring(0, host.length()-1) : host
        }
        catch( URISyntaxException e ) {
            return null
        }
    }

    /**
     * Whether a {@code provider/model} identifier targets the OpenAI wire protocol. Protocol only:
     * the credential and the endpoint follow {@link #apiProviderFor} for every provider.
     */
    static boolean isOpenAiProtocol(String modelId) {
        return providerPrefixOf(modelId) == OPENAI_PROVIDER
    }

    /** The lower-case provider prefix of a {@code provider/model} id, or {@code null}. */
    static String providerPrefixOf(String modelId) {
        if( !modelId )
            return null
        final i = modelId.indexOf('/')
        return i > 0 ? modelId.substring(0, i).toLowerCase() : null
    }

    /**
     * The provider whose namespace the credential and the endpoint are resolved from for the given
     * {@code provider/model} id (design D1): the explicit {@code agent.apiProvider}, else the
     * provider inferred from the provider-neutral endpoint, else the model-id prefix -- today's
     * answer, and the right one for a runner (pi) whose prefix already names a real provider.
     *
     * <p>Inference outranks the prefix because the prefix names a PROTOCOL: {@code openai/} with a
     * {@code baseUrl} of {@code https://openrouter.ai/api/v1} is the documented way to reach
     * OpenRouter, and the credential that endpoint wants is OpenRouter's.
     *
     * <p>A method rather than a field because the effective model is the per-agent {@code model}
     * directive joined with the {@code agent.model} default one layer up (see
     * {@link nextflow.script.AgentDef#buildAgentTask}) -- this class never sees it. The two rungs
     * that do not need it ARE fields.
     */
    String apiProviderFor(String modelId) {
        return apiProvider ?: inferredProvider ?: providerPrefixOf(modelId)
    }

    /**
     * The credential to hand to a runner for the given {@code provider/model} id: the neutral
     * tiers, else the provider's own variable WHEN the endpoint it would travel to belongs to that
     * provider.
     *
     * <p>The neutral tiers carry no assumption -- the user named that value for THIS agent,
     * whatever it targets -- so they always win.
     *
     * <p>The provider tier is different: it is a variable exported for a PROVIDER, not for this
     * endpoint, and a runner installs whatever it is handed as the credential OF THE MODEL'S
     * PROVIDER, ahead of everything it could resolve itself (pi's {@code setRuntimeApiKey} override
     * is read before its credential store and before the ambient environment). Sending it to an
     * arbitrary address is therefore a credential disclosure the pipeline chooses, so it travels
     * only when the endpoint agrees: no endpoint at all AND a provider that is the model prefix's
     * own (so the default endpoint dialled is that provider's), the provider's own host, the
     * {@code <PROVIDER>_BASE_URL} that same namespace supplied, or an endpoint vouched for by an
     * explicit {@code agent.apiProvider}. A corporate gateway or another provider's host gets
     * nothing -- and is told what to set.
     */
    String apiKeyFor(String modelId) {
        if( apiKey )
            return apiKey
        final provider = apiProviderFor(modelId)
        // D3 observability: the prefix is the answer a reader assumes, so name the other one when
        // it wins -- otherwise which variable gets read is invisible
        final prefix = providerPrefixOf(modelId)
        if( provider && prefix && provider != prefix )
            log.debug "Resolved API provider `${provider}` for agent model `${modelId}` from ${apiProvider ? '`agent.apiProvider`' : 'the endpoint ' + baseUrl} (the model prefix is `${prefix}`)"
        final resolved = resolveProviderApiKey(sysEnv, provider)
        if( !resolved )
            return null
        if( !providerTierWithheld(modelId) )
            return resolved
        final endpoint = baseUrlFor(modelId)
        log.warn withheldCredentialWarning(modelId, provider, prefix, endpoint)
        return null
    }

    /**
     * The WARN {@link #apiKeyFor} emits when the provider tier resolved a credential and the
     * endpoint gate then withheld it: which variable was read, why its value is not being sent,
     * and the two ways to proceed.
     *
     * <p>Not {@code static}: the variable name is the one {@code sysEnv} actually set (a provider
     * with an alias pair can be read through either), so the message cannot be derived from
     * {@code provider} alone.
     *
     * @param modelId  the {@code provider/model} id the credential was resolved for
     * @param provider the namespace the credential came from, per {@link #apiProviderFor}
     * @param prefix   the model id's own {@code provider/} prefix, or null when it has none
     * @param endpoint the resolved endpoint, or null when none resolved
     */
    private String withheldCredentialWarning(String modelId, String provider, String prefix, String endpoint) {
        final String because
        if( endpoint )
            because = "the resolved endpoint ${endpoint} is not a known `${provider}` endpoint"
        else if( prefix )
            because = "no endpoint resolved, so the request goes to the default endpoint of the model's own `${prefix}` provider and not to `${provider}`"
        else
            // a model id with no `provider/` prefix names no default endpoint either, so there is
            // nothing this credential can be shown to belong to
            because = "no endpoint resolved and the model id names no provider, so there is no endpoint this credential is known to belong to"
        final remedy = endpoint
            ? "`agent.apiProvider = '${provider}'` to confirm it accepts `${provider}` credentials"
            : "`agent.baseUrl` to the `${provider}` endpoint this credential belongs to"
        return "Not using the ${providerApiKeyVar(sysEnv, provider)} credential for agent model `${modelId}`: ${because} - set `agent.apiKey` to the credential that endpoint accepts, or ${remedy}"
    }

    /**
     * Whether a {@code <PROVIDER>_API_KEY} DID resolve for this model and was then withheld by the
     * endpoint gate, as opposed to nothing resolving at all. The two look identical through
     * {@link #apiKeyFor}, which answers {@code null} to both, yet they are opposite situations: the
     * second is the ordinary "no credential here" that a local endpoint (or a runner with sources
     * of its own) handles perfectly well, while the first is a MISCONFIGURATION whose credential is
     * sitting right there, unusable. Only the first must be reported as an error by a runner that
     * has no other source, and neither may be turned into
     * {@link AgentRunnerRequest#PLACEHOLDER_API_KEY} -- a placeholder buys an opaque 401 in place of
     * a diagnosis that is available.
     *
     * <p>Silent by design: {@link #apiKeyFor} already emits the WARN, and this is queried beside it
     * for the same model.
     */
    boolean credentialWithheldFor(String modelId) {
        return providerTierWithheld(modelId)
    }

    /**
     * The gate decision, without the logging: {@code true} when the provider tier resolved a value
     * that {@link #providerCredentialApplies} refuses to send. Shared by {@link #apiKeyFor} and
     * {@link #credentialWithheldFor} so the credential and the signal describing it can never
     * disagree.
     */
    private boolean providerTierWithheld(String modelId) {
        if( apiKey )
            return false
        final provider = apiProviderFor(modelId)
        if( !resolveProviderApiKey(sysEnv, provider) )
            return false
        return !providerCredentialApplies(provider, baseUrlFor(modelId), modelId)
    }

    /**
     * Whether the {@code <PROVIDER>_API_KEY} tier may travel to the endpoint resolved for this
     * agent -- see {@link #apiKeyFor} for why an ambient provider variable needs a matching
     * endpoint and the neutral tiers do not.
     */
    private boolean providerCredentialApplies(String provider, String endpoint, String modelId) {
        if( !endpoint )
            // No endpoint resolved does NOT mean "the resolved provider's own default". The runner
            // dials the default endpoint of the provider named by the MODEL PREFIX -- the
            // langchain4j client hardcodes https://api.openai.com/v1, and pi dials the provider its
            // prefix names -- so a credential resolved in a DIFFERENT namespace would be sent to a
            // third party. Concretely: `agent.apiProvider = 'openrouter'` with `openai/gpt-4o` and
            // no `agent.baseUrl` used to ship OPENROUTER_API_KEY to OpenAI.
            return provider == providerPrefixOf(modelId)
        final endpointProvider = inferProviderFromUrl(endpoint)
        if( endpointProvider != null )
            // a host this table knows: only its own provider's credential belongs there, whatever
            // the config says -- this is also what stops OPENAI_BASE_URL=https://openrouter.ai/...
            // from resolving the openai credential by way of the model prefix
            return endpointProvider == provider
        // an unrecognized host (a gateway, a local server): the variable was not exported for it,
        // so it travels only with a statement that it may
        return apiProvider != null || endpoint == resolveProviderBaseUrl(sysEnv, provider)
    }

    /**
     * The endpoint to hand to a runner for the given {@code provider/model} id: the neutral tiers,
     * else the provider's own {@code <PROVIDER>_BASE_URL}. NOT gated like {@link #apiKeyFor} -- an
     * endpoint is not a secret, and a variable named for a provider, used for that provider's
     * model, is the intent rather than an accident.
     *
     * <p>NOTE: the resolved value enters the resume cache key (see
     * {@link nextflow.script.AgentDef#canonicalAgentSource}), so a run that already exports one of
     * these variables sees a one-time cache invalidation for the agents it now resolves for.
     */
    String baseUrlFor(String modelId) {
        return baseUrl ?: resolveProviderBaseUrl(sysEnv, apiProviderFor(modelId))
    }

    /**
     * The remedy to name when no credential resolved for a provider: the variables the ladder
     * ACTUALLY consults for it. A message naming {@code OPENAI_API_KEY} while the run resolved an
     * {@code anthropic} provider is worse than no message.
     */
    static String missingCredentialHint(String provider) {
        final names = new ArrayList<String>()
        names.add('NXF_AGENT_API_KEY')
        names.addAll(apiKeyVarsFor(provider))
        final vars = names.collect { "`${it}`" }.join(' or ')
        final hint = "set `agent.apiKey` in the Nextflow configuration, or the ${vars} environment variable"
        return isKnownProvider(provider) ? hint : hint + ", and name the credential namespace with `agent.apiProvider`"
    }

    String getRunner() { runner }

    String getModel() { model }

    String getApiProvider() { apiProvider }

    /** The PROVIDER-NEUTRAL credential; what a runner gets is {@link #apiKeyFor}. */
    String getApiKey() { apiKey }

    /** The PROVIDER-NEUTRAL endpoint; what a runner gets is {@link #baseUrlFor}. */
    String getBaseUrl() { baseUrl }

    Integer getMaxIterations() { maxIterations }

    Duration getRequestTimeout() { requestTimeout }

    MemoryUnit getMaxToolOutputInlineSize() { maxToolOutputInlineSize }

    Boolean getTrace() { trace }

    AgentRpcConfig getRpc() { rpc }

    /**
     * The effective maximum size, in bytes, of a structured tool-output file whose contents
     * are inlined for the LLM; defaults to 32 KB when not configured.
     */
    long maxToolOutputInlineBytes() { maxToolOutputInlineSize != null ? maxToolOutputInlineSize.toBytes() : ToolOutputReader.DEFAULT_INLINE_BYTES }

    /** Whether execution tracing is enabled (defaults to {@code false} when not configured). */
    boolean traceEnabled() { trace != null && trace }
}
