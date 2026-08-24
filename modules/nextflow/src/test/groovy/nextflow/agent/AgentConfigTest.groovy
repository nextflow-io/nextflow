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

import nextflow.SysEnv
import nextflow.config.ConfigValidator
import nextflow.exception.AbortOperationException
import nextflow.script.dsl.ProcessBuilder
import nextflow.util.Duration
import org.junit.Rule
import spock.lang.Specification
import test.OutputCapture
import nextflow.agent.rpc.AgentRpcHostResolver
import nextflow.agent.rpc.AgentRpcConfig

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AgentConfigTest extends Specification {

    @Rule
    OutputCapture capture = new OutputCapture()

    def 'should build from a config map'() {
        when:
        def config = new AgentConfig([runner: 'pi', model: 'openai/gpt-5-mini', maxIterations: 7, requestTimeout: '90s', maxToolOutputInlineSize: '64 KB', rpc: [port: 1234, remoteHost: 'driver.internal']])
        then:
        config.runner == 'pi'
        config.model == 'openai/gpt-5-mini'
        config.maxIterations == 7
        config.requestTimeout == Duration.of('90s')
        config.maxToolOutputInlineBytes() == 64 * 1024
        config.rpc.port == 1234
        config.rpc.remoteHost == 'driver.internal'
    }

    def 'should default the RPC capability timeout to one hour'() {
        expect: 'the default absorbs an executor queueing delay -- the clock starts when the task SCRIPT is generated'
        AgentRpcConfig.DEFAULT_CAPABILITY_TIMEOUT == Duration.of('1h')
        new AgentConfig([:]).rpc.capabilityTimeout == Duration.of('1h')

        and: 'an operator can widen or tighten it, from a string or a Duration'
        new AgentConfig([rpc: [capabilityTimeout: '10m']]).rpc.capabilityTimeout == Duration.of('10m')
        new AgentConfig([rpc: [capabilityTimeout: Duration.of('2h')]]).rpc.capabilityTimeout == Duration.of('2h')
    }

    def 'should accept agent.rpc.capabilityTimeout as a declared config option'() {
        when:
        new ConfigValidator().validate([agent: [rpc: [capabilityTimeout: '10m']]])
        then:
        !capture.toString().contains("Unrecognized config option 'agent")
    }

    def 'should enable RPC transport security unless it is explicitly disabled'() {
        expect: 'secure by default, including for the no-arg form used by the spec reflection'
        new AgentConfig([:]).rpc.tlsEnabled()
        new AgentRpcConfig().tlsEnabled()

        and: 'only an explicit false opts out'
        !new AgentConfig([rpc: [tls: false]]).rpc.tlsEnabled()
        new AgentConfig([rpc: [tls: true]]).rpc.tlsEnabled()

        when:
        new ConfigValidator().validate([agent: [rpc: [tls: false]]])
        then:
        !capture.toString().contains("Unrecognized config option 'agent")
    }

    def 'should build from an empty config map'() {
        when:
        def config = new AgentConfig([:])
        then:
        AgentConfig.DEFAULT_EXECUTOR == 'local'
        config.model == null
        config.maxIterations == null
        config.requestTimeout == null
        config.maxToolOutputInlineBytes() == 32768
        and:
        config.trace == null
        !config.traceEnabled()
        config.rpc.port == AgentRpcConfig.DEFAULT_PORT
        and: 'the broker host is NOT defaulted: only the container engine can stand in for it'
        config.rpc.remoteHost == null
    }

    def 'the broker host is inferred from the engine when it is not configured'() {
        given:
        def unset = new AgentRpcConfig([:])
        def explicit = new AgentRpcConfig([remoteHost: 'driver.internal'])

        expect: 'docker and podman name the container host themselves'
        unset.resolveRemoteHost('docker') == 'host.docker.internal'
        unset.resolveRemoteHost('podman') == 'host.containers.internal'

        and: 'an engine that creates no network namespace shares the driver`s, so loopback reaches it'
        unset.resolveRemoteHost('singularity') == '127.0.0.1'

        and: 'no engine at all is still nothing to go on -- see AgentRpcHostResolver error row E7'
        unset.resolveRemoteHost(null) == null

        and: 'an explicit value always wins over every inferred row'
        explicit.resolveRemoteHost('docker') == 'driver.internal'
        explicit.resolveRemoteHost(null) == 'driver.internal'
    }

    def 'the broker host resolves config, then the environment, then the engine alias'() {
        expect: 'each rung is reached only when the ones above it are absent'
        AgentRpcConfig.resolveConfiguredHost(opts, env) == expected

        where:
        opts                            | env                                                    || expected
        [remoteHost: 'from.config']     | [:]                                                    || 'from.config'
        [:]                             | [NXF_AGENT_RPC_REMOTE_HOST: 'from.env']                || 'from.env'
        // config wins: it is written for THIS pipeline, the variable for whatever shares the shell
        [remoteHost: 'from.config']     | [NXF_AGENT_RPC_REMOTE_HOST: 'from.env']                || 'from.config'
        // neither set leaves the engine alias to answer, via resolveRemoteHost
        [:]                             | [:]                                                    || null
        // an empty value falls THROUGH rather than shadowing the rung below and advertising ''
        [remoteHost: '']                | [NXF_AGENT_RPC_REMOTE_HOST: 'from.env']                || 'from.env'
        [:]                             | [NXF_AGENT_RPC_REMOTE_HOST: '']                        || null
        [remoteHost: '']                | [NXF_AGENT_RPC_REMOTE_HOST: '']                        || null
    }

    def 'the environment reaches the config object, and the engine alias still backstops it'() {
        given:
        SysEnv.push([NXF_AGENT_RPC_REMOTE_HOST: 'driver.from.env'])

        when: 'nothing is configured'
        def fromEnv = new AgentRpcConfig([:])

        then: 'the variable is what both AgentDef and the broker read'
        fromEnv.remoteHost == 'driver.from.env'
        and: 'and it outranks the engine alias, being explicit about the driver'
        fromEnv.resolveRemoteHost('docker') == 'driver.from.env'

        when: 'the pipeline configures its own'
        def fromConfig = new AgentRpcConfig([remoteHost: 'driver.from.config'])

        then:
        fromConfig.remoteHost == 'driver.from.config'

        cleanup:
        SysEnv.pop()
    }

    def 'only docker and podman name the host their containers run on'() {
        expect: 'the alias table, which is ONE row of the ladder AgentRpcHostResolver owns'
        AgentRpcConfig.hostAliasFor('docker') == 'host.docker.internal'
        AgentRpcConfig.hostAliasFor('podman') == 'host.containers.internal'

        and: 'no other engine has one -- which is now a reason to infer, not a reason to reject'
        for( final engine : ['singularity', 'apptainer', 'sarus', 'shifter', 'charliecloud', 'apple-container', 'smolvm'] )
            assert AgentRpcConfig.hostAliasFor(engine) == null

        and: 'an unknown or absent engine is not an error, just no alias'
        AgentRpcConfig.hostAliasFor('nonesuch') == null
        AgentRpcConfig.hostAliasFor(null) == null
    }

    def 'should read the trace option'() {
        expect:
        new AgentConfig([trace: true]).traceEnabled()
        and:
        !new AgentConfig([trace: false]).traceEnabled()
    }

    def 'AGENT_ONLY_OPTIONS is derived from the declared options and disjoint from process directives'() {
        expect: 'the agent-only axis is the @ConfigOption fields PLUS the nested config scopes'
        AgentConfig.AGENT_ONLY_OPTIONS == [
            'runner', 'model', 'apiProvider', 'apiKey', 'baseUrl', 'maxIterations', 'requestTimeout',
            'maxToolOutputInlineSize', 'trace', 'rpc' ] as Set

        and: 'drift guard: re-adding a task directive as an @ConfigOption would break the two axes apart'
        AgentConfig.AGENT_ONLY_OPTIONS.intersect(ProcessBuilder.DIRECTIVES as Set).isEmpty()
    }

    def 'resolveOptions applies the plain scope only'() {
        when:
        def opts = AgentConfig.resolveOptions([model: 'm0', cpus: 4], [], 'qa', 'qa', 'qa')
        then: 'agent-only options are copied; task directives are not'
        opts == [model: 'm0']
    }

    def 'resolveOptions applies the selector ladder, weakest to strongest'() {
        given:
        def scope = [
            model: 'plain',
            'withLabel:big': [model: 'label'],
            'withName:qa': [model: 'base'],
            'withName:reviewer': [model: 'alias'],
            'withName:WF:qa': [model: 'fq'],
        ]

        and:
        def without = { String... keys -> scope.findAll { e -> !(e.key in (keys as List)) } }

        expect: 'fully-qualified wins over alias, base, label and plain'
        AgentConfig.resolveOptions(scope, ['big'], 'qa', 'reviewer', 'WF:qa').model == 'fq'
        and: 'without the fq rule the alias wins'
        AgentConfig.resolveOptions(without('withName:WF:qa'), ['big'], 'qa', 'reviewer', 'WF:qa').model == 'alias'
        and: 'without the alias rule the base name wins'
        AgentConfig.resolveOptions(without('withName:WF:qa', 'withName:reviewer'), ['big'], 'qa', 'reviewer', 'WF:qa').model == 'base'
        and: 'without any name rule the label rule wins'
        AgentConfig.resolveOptions([model: 'plain', 'withLabel:big': [model: 'label']], ['big'], 'qa', 'qa', 'qa').model == 'label'
        and: 'with no matching selector the plain scope applies'
        AgentConfig.resolveOptions(scope, [], 'other', 'other', 'other').model == 'plain'
    }

    def 'resolveOptions matches regex and negated selector targets'() {
        expect: 'a regex withName target'
        AgentConfig.resolveOptions(['withName:shard.+': [model: 'm']], [], 'shard_1', 'shard_1', 'shard_1').model == 'm'
        and: 'a negated withName target does not match the excluded agent'
        AgentConfig.resolveOptions(['withName:!critic': [model: 'm']], [], 'critic', 'critic', 'critic').isEmpty()
        and: 'but does match any other agent'
        AgentConfig.resolveOptions(['withName:!critic': [model: 'm']], [], 'qa', 'qa', 'qa').model == 'm'
        and: 'a negated withLabel target matches a label-less agent'
        AgentConfig.resolveOptions(['withLabel:!big': [model: 'm']], [], 'qa', 'qa', 'qa').model == 'm'
    }

    def 'resolveOptions ignores task directives and tolerates a malformed selector body'() {
        expect: 'a task directive inside a selector is NOT copied onto the agent-only axis'
        AgentConfig.resolveOptions(['withName:qa': [cpus: 8, model: 'm']], [], 'qa', 'qa', 'qa') == [model: 'm']

        and: 'a selector-scoped nested scope IS copied (a documented no-op: the broker reads the session config)'
        AgentConfig.resolveOptions(['withName:qa': [rpc: [port: 1234]]], [], 'qa', 'qa', 'qa') == [rpc: [port: 1234]]

        and: 'a non-map selector body is skipped without a ClassCastException (applyConfig reports it)'
        AgentConfig.resolveOptions(['withName:qa': 'oops'], [], 'qa', 'qa', 'qa').isEmpty()
    }

    def 'should be a recognized config scope'() {
        when:
        new ConfigValidator().validate([
            agent: [
                runner: 'pi',
                executor: 'k8s',
                container: 'runner:1',
                arch: 'arm64',
                cpus: 2,
                memory: '1 GB',
                time: '1h',
                queue: 'agents',
                errorStrategy: 'retry',
                maxRetries: 2,
                cache: false,
                model: 'openai/gpt-5-mini',
                apiProvider: 'openai',
                apiKey: 'sk-xxx',
                baseUrl: 'http://localhost:8000/v1',
                maxIterations: 7,
                requestTimeout: '90s',
                maxToolOutputInlineSize: '64 KB',
                trace: true,
                rpc: [port: 1234, remoteHost: 'driver.internal']
            ]
        ])
        then:
        !capture.toString().contains("Unrecognized config option 'agent")
    }

    // -----------------------------------------------------------------------
    // Endpoint and credential resolution (design D1/D2/D3): ONE ladder for every provider --
    // config, then the provider-neutral NXF_AGENT_*, then the PROVIDER's own variable -- resolved
    // here in core so neither runner reads the environment. Driven with a FAKE env map (the
    // two-arg constructor, the house pattern of `AwsConfig.getAwsRegion`) so the real process
    // environment is never mutated.
    // -----------------------------------------------------------------------

    def 'resolveNeutralApiKey applies the config -> NXF_AGENT_API_KEY tiers'() {
        expect: 'nothing set resolves to null: whether a credential is REQUIRED is a runner concern'
        AgentConfig.resolveNeutralApiKey([:], [:]) == null

        and: 'tier 1 is the config option (already selector-merged into opts by resolveOptions)'
        AgentConfig.resolveNeutralApiKey([:], [apiKey: 'sk-config']) == 'sk-config'

        and: 'the config option beats the environment'
        AgentConfig.resolveNeutralApiKey([NXF_AGENT_API_KEY: 'sk-nxf'], [apiKey: 'sk-config']) == 'sk-config'

        and: 'tier 2 is the provider-neutral variable'
        AgentConfig.resolveNeutralApiKey([NXF_AGENT_API_KEY: 'sk-nxf'], [:]) == 'sk-nxf'

        and: 'the PROVIDER tier is not reachable from here: it needs the model, which this does not see'
        AgentConfig.resolveNeutralApiKey([OPENAI_API_KEY: 'sk-openai'], [:]) == null

        and: 'a null env or a null opts map is tolerated (the no-arg extension-point ctor)'
        AgentConfig.resolveNeutralApiKey(null, null) == null
        AgentConfig.resolveNeutralApiKey(null, [apiKey: 'sk-config']) == 'sk-config'
        AgentConfig.resolveNeutralApiKey([NXF_AGENT_API_KEY: 'sk-nxf'], null) == 'sk-nxf'
    }

    def 'resolveNeutralBaseUrl applies the config -> NXF_AGENT_BASE_URL tiers'() {
        expect: 'nothing set resolves to null, which means "use the provider default"'
        AgentConfig.resolveNeutralBaseUrl([:], [:]) == null

        and: 'tier 1 is the config option, and it beats the environment'
        AgentConfig.resolveNeutralBaseUrl([:], [baseUrl: 'http://config:8000/v1']) == 'http://config:8000/v1'
        AgentConfig.resolveNeutralBaseUrl([NXF_AGENT_BASE_URL: 'http://nxf/v1'], [baseUrl: 'http://config:8000/v1']) == 'http://config:8000/v1'

        and: 'tier 2 is the provider-neutral variable'
        AgentConfig.resolveNeutralBaseUrl([NXF_AGENT_BASE_URL: 'http://nxf/v1'], [:]) == 'http://nxf/v1'

        and: 'and the provider tier is deliberately absent here -- it is what the inference feeds on (D3)'
        AgentConfig.resolveNeutralBaseUrl([OPENAI_BASE_URL: 'http://openai/v1'], [:]) == null

        and: 'a null env or a null opts map is tolerated'
        AgentConfig.resolveNeutralBaseUrl(null, null) == null
    }

    def 'the ladder is truthiness-based so an empty value does not shadow the tiers below it'() {
        given: 'an EMPTY value is not a credential: `agent.apiKey = params.key` with `params.key`'
        // unset, or `export NXF_AGENT_API_KEY=`, both yield ''. A null-check ladder would let that
        // win its tier, shadow everything below and then trip the runner's "missing credential"
        // branch with a misleading message. Consequence worth knowing: an explicitly empty
        // `agent.apiKey` does NOT mean "no credential" -- the D8 no-credential path is reached by
        // leaving the option unset while setting a baseUrl, not by setting it empty.
        expect:
        AgentConfig.resolveNeutralApiKey([NXF_AGENT_API_KEY: 'sk-nxf'], [apiKey: '']) == 'sk-nxf'
        AgentConfig.resolveNeutralBaseUrl([NXF_AGENT_BASE_URL: 'http://nxf/v1'], [baseUrl: '']) == 'http://nxf/v1'

        and: 'with nothing to fall back on an empty value normalizes to null, not to ""'
        AgentConfig.resolveNeutralApiKey([:], [apiKey: '']) == null
        AgentConfig.resolveNeutralBaseUrl([:], [baseUrl: '']) == null

        and: 'the PROVIDER tier applies the same rule: an exported-but-empty variable falls through'
        AgentConfig.resolveProviderApiKey([GEMINI_API_KEY: '', GOOGLE_API_KEY: 'sk-goo'], 'gemini') == 'sk-goo'
        AgentConfig.resolveProviderApiKey([OPENAI_API_KEY: ''], 'openai') == null
        AgentConfig.resolveProviderBaseUrl([OPENAI_BASE_URL: ''], 'openai') == null

        and: 'and an empty value reaches the provider tier through the whole ladder'
        new AgentConfig([apiKey: ''], [NXF_AGENT_API_KEY: '', OPENAI_API_KEY: 'sk-openai'])
            .apiKeyFor('openai/gpt-5-mini') == 'sk-openai'
    }

    def 'the constructor resolves the NEUTRAL tiers through SysEnv'() {
        given: 'resolution happens in the CONSTRUCTOR (as in AwsConfig), reading SysEnv and never System.getenv'
        SysEnv.push(ENV)
        and:
        def config = new AgentConfig(OPTS)

        expect: 'the FIELDS carry the provider-neutral tiers only -- what `nextflow config` shows'
        // The provider tier is per-MODEL (which variable is read depends on the effective model),
        // so it cannot be a field; folding it in here would also make `nextflow config` MORE
        // revealing about the ambient environment than it is today.
        config.apiKey == API_KEY
        config.baseUrl == BASE_URL

        cleanup:
        SysEnv.pop()

        where:
        ENV                                                                | OPTS                                             || API_KEY     | BASE_URL
        [:]                                                                | [:]                                              || null        | null
        [:]                                                                | [apiKey: 'sk-config', baseUrl: 'http://cfg/v1']  || 'sk-config' | 'http://cfg/v1'
        [NXF_AGENT_API_KEY: 'sk-nxf', NXF_AGENT_BASE_URL: 'http://nxf/v1'] | [:]                                              || 'sk-nxf'    | 'http://nxf/v1'
        [NXF_AGENT_API_KEY: 'sk-nxf', OPENAI_API_KEY: 'sk-oai']            | [:]                                              || 'sk-nxf'    | null
        [OPENAI_API_KEY: 'sk-oai']                                         | [apiKey: 'sk-config']                            || 'sk-config' | null
        and: 'the provider tier is NOT a field: it is resolved per model by apiKeyFor/baseUrlFor'
        [OPENAI_API_KEY: 'sk-oai', OPENAI_BASE_URL: 'http://oai/v1']       | [:]                                              || null        | null
        and: 'the two options resolve INDEPENDENTLY: a local endpoint from config, the key from the env'
        [NXF_AGENT_API_KEY: 'sk-nxf']                                      | [baseUrl: 'http://localhost:8000/v1']            || 'sk-nxf'    | 'http://localhost:8000/v1'
        and: 'a local endpoint and NO credential at all is a valid combination (design D8)'
        [:]                                                                | [baseUrl: 'http://localhost:11434/v1']           || null        | 'http://localhost:11434/v1'
    }

    def 'the two-arg constructor pins the environment the deferred provider tier reads'() {
        given: 'the provider tier resolves per MODEL, long after construction, so the environment is'
        // HELD rather than re-read: an object whose answers depend on when they are asked cannot be
        // reasoned about. The explicit form is also what lets a spec swap a map instead of SysEnv.
        SysEnv.push([OPENAI_API_KEY: 'sk-from-sysenv'])
        and:
        def config = new AgentConfig([:], [OPENAI_API_KEY: 'sk-injected'])

        expect:
        config.apiKeyFor('openai/gpt-5-mini') == 'sk-injected'

        cleanup:
        SysEnv.pop()
    }

    def 'the provider-neutral tiers reach every provider and are never endpoint-gated'() {
        given: 'the config option and NXF_AGENT_* are named by the user FOR THIS AGENT, whichever'
        // provider it targets and whatever it points at, so unlike an ambient `<PROVIDER>_API_KEY`
        // they carry no assumption that has to be checked against the endpoint
        def env = [NXF_AGENT_API_KEY: 'sk-nxf', NXF_AGENT_BASE_URL: 'http://nxf/v1', OPENAI_API_KEY: 'sk-openai']

        expect: 'tier 2 applies to a non-openai model'
        with(new AgentConfig([:], env)) {
            apiKeyFor('anthropic/claude-sonnet-4') == 'sk-nxf'
            baseUrlFor('anthropic/claude-sonnet-4') == 'http://nxf/v1'
        }

        and: 'and so does tier 1, which also beats it'
        with(new AgentConfig([apiKey: 'sk-config', baseUrl: 'http://config/v1'], env)) {
            apiKeyFor('anthropic/claude-sonnet-4') == 'sk-config'
            baseUrlFor('anthropic/claude-sonnet-4') == 'http://config/v1'
            apiKeyFor('openai/gpt-5-mini') == 'sk-config'
        }

        and: 'a neutral credential travels to ANOTHER provider\'s own host without a murmur'
        new AgentConfig([apiKey: 'sk-config', baseUrl: 'https://api.anthropic.com/v1'], [:])
            .apiKeyFor('openai/gpt-5-mini') == 'sk-config'
    }

    // -----------------------------------------------------------------------
    // D1: which provider namespace the credential and the endpoint come from. Three rungs, and
    // NONE of them selects the wire protocol -- that stays the model-id prefix.
    // -----------------------------------------------------------------------

    def 'apiProviderFor resolves explicit, then the inferred endpoint, then the model prefix'() {
        expect: 'the prefix is the historical answer, and still the last word'
        new AgentConfig([:], [:]).apiProviderFor('anthropic/claude-sonnet-4') == 'anthropic'

        and: 'a provider-neutral endpoint on a well-known host OUTRANKS it -- the prefix names a'
        // PROTOCOL, so `openai/` + https://openrouter.ai/api/v1 is the documented way to reach
        // OpenRouter and the credential that endpoint wants is OpenRouter's
        new AgentConfig([baseUrl: 'https://openrouter.ai/api/v1'], [:]).apiProviderFor('openai/gpt-4o') == 'openrouter'
        new AgentConfig([:], [NXF_AGENT_BASE_URL: 'https://api.mistral.ai/v1']).apiProviderFor('openai/x') == 'mistral'

        and: 'an explicit `agent.apiProvider` outranks both'
        new AgentConfig([apiProvider: 'azure', baseUrl: 'https://openrouter.ai/api/v1'], [:])
            .apiProviderFor('openai/gpt-4o') == 'azure'

        and: 'an unrecognized endpoint host does not fire, so the prefix answers'
        new AgentConfig([baseUrl: 'https://gateway.corp/v1'], [:]).apiProviderFor('openai/gpt-4o') == 'openai'

        and: 'inference reads the NEUTRAL endpoint only: `<PROVIDER>_BASE_URL` cannot feed the'
        // inference that decides which provider's variable to read in the first place (D3 circularity)
        new AgentConfig([:], [OPENAI_BASE_URL: 'https://openrouter.ai/api/v1']).apiProviderFor('openai/x') == 'openai'

        and: 'with no prefix and nothing to infer from there is no provider at all'
        new AgentConfig([:], [:]).apiProviderFor('gpt-5-mini') == null
        new AgentConfig([:], [:]).apiProviderFor(null) == null
    }

    def 'providerPrefixOf reads the model-id prefix and isOpenAiProtocol is exactly that prefix'() {
        expect:
        AgentConfig.providerPrefixOf('openai/gpt-5-mini') == 'openai'
        AgentConfig.providerPrefixOf('Anthropic/Claude-Sonnet-4') == 'anthropic'
        AgentConfig.providerPrefixOf('openrouter/openai/gpt-4o') == 'openrouter'
        and: 'no prefix is not an error, just no answer'
        AgentConfig.providerPrefixOf('gpt-5-mini') == null
        AgentConfig.providerPrefixOf('/gpt-5-mini') == null
        AgentConfig.providerPrefixOf('') == null
        AgentConfig.providerPrefixOf(null) == null

        and: 'the protocol predicate is the prefix and nothing else -- it no longer gates the ladder'
        AgentConfig.isOpenAiProtocol('openai/gpt-5-mini')
        AgentConfig.isOpenAiProtocol('openai/llama-3.3-70b')
        !AgentConfig.isOpenAiProtocol('anthropic/claude-sonnet-4')
        !AgentConfig.isOpenAiProtocol('openrouter/openai/gpt-4o')
        !AgentConfig.isOpenAiProtocol('gpt-5-mini')
        !AgentConfig.isOpenAiProtocol(null)
    }

    def 'agent.apiProvider is normalized, and an unknown value is rejected rather than uppercased'() {
        expect: 'unset (or empty, by the same truthiness rule as every other tier) is null'
        AgentConfig.resolveApiProvider([:]) == null
        AgentConfig.resolveApiProvider([apiProvider: '']) == null
        AgentConfig.resolveApiProvider([apiProvider: '   ']) == null
        AgentConfig.resolveApiProvider(null) == null

        and: 'trimmed and lower-cased, but never otherwise mangled'
        AgentConfig.resolveApiProvider([apiProvider: ' OpenAI ']) == 'openai'
        AgentConfig.resolveApiProvider([apiProvider: 'ANTHROPIC']) == 'anthropic'

        when: 'an unrecognized token is written'
        new AgentConfig([apiProvider: 'stripe'], [:])

        then: 'it fails with the accepted values, rather than silently reading no variable at all'
        // the closed namespace is what keeps `<PROVIDER>_API_KEY` from naming an arbitrary variable
        // in the driver's environment; an unrecognized token would otherwise be a typo diagnosed nowhere
        def error = thrown(AbortOperationException)
        error.message.contains('`agent.apiProvider`')
        error.message.contains('stripe')
        error.message.contains('anthropic, azure, gemini, google, mistral, openai, openrouter')
    }

    def 'the provider namespace is a closed, explicit table'() {
        expect: 'every token the ladder knows, and the variables it reads for each, in order'
        // NOTE: the key names must stay equal to `AgentSecretMasker.SECRET_ENV_KEYS` minus the
        // neutral NXF_AGENT_API_KEY -- the redaction backstop and the resolution contract must not
        // disagree about what counts as a credential.
        AgentConfig.knownProviders() == ['anthropic', 'azure', 'gemini', 'google', 'mistral', 'openai', 'openrouter'] as Set
        AgentConfig.apiKeyVarsFor('anthropic') == ['ANTHROPIC_API_KEY']
        AgentConfig.apiKeyVarsFor('azure') == ['AZURE_OPENAI_API_KEY']
        AgentConfig.apiKeyVarsFor('gemini') == ['GEMINI_API_KEY', 'GOOGLE_API_KEY']
        AgentConfig.apiKeyVarsFor('google') == ['GOOGLE_API_KEY', 'GEMINI_API_KEY']
        AgentConfig.apiKeyVarsFor('mistral') == ['MISTRAL_API_KEY']
        AgentConfig.apiKeyVarsFor('openai') == ['OPENAI_API_KEY']
        AgentConfig.apiKeyVarsFor('openrouter') == ['OPENROUTER_API_KEY']

        and: 'the endpoint half is NOT the same key set: only these three variables actually exist'
        // reading a MISTRAL_BASE_URL would be Nextflow inventing a convention under a vendor's name
        AgentConfig.baseUrlVarsFor('openai') == ['OPENAI_BASE_URL']
        AgentConfig.baseUrlVarsFor('anthropic') == ['ANTHROPIC_BASE_URL']
        AgentConfig.baseUrlVarsFor('azure') == ['AZURE_OPENAI_ENDPOINT']
        AgentConfig.baseUrlVarsFor('gemini') == []
        AgentConfig.baseUrlVarsFor('google') == []
        AgentConfig.baseUrlVarsFor('mistral') == []
        AgentConfig.baseUrlVarsFor('openrouter') == []

        and: 'an unknown token names NO variable, so a model prefix can never reach one'
        !AgentConfig.isKnownProvider('stripe')
        !AgentConfig.isKnownProvider(null)
        AgentConfig.apiKeyVarsFor('stripe') == []
        AgentConfig.baseUrlVarsFor('stripe') == []

        and: 'the candidate list is ordered: the first TRUTHY hit wins, not the last'
        AgentConfig.providerApiKeyVar([GEMINI_API_KEY: 'a', GOOGLE_API_KEY: 'b'], 'gemini') == 'GEMINI_API_KEY'
        AgentConfig.providerApiKeyVar([GEMINI_API_KEY: 'a', GOOGLE_API_KEY: 'b'], 'google') == 'GOOGLE_API_KEY'
        and: 'and the alias answers when the preferred spelling is absent'
        AgentConfig.providerApiKeyVar([GOOGLE_API_KEY: 'b'], 'gemini') == 'GOOGLE_API_KEY'
        AgentConfig.providerApiKeyVar([GEMINI_API_KEY: 'a'], 'google') == 'GEMINI_API_KEY'
        AgentConfig.providerApiKeyVar([:], 'gemini') == null
        AgentConfig.providerApiKeyVar(null, 'gemini') == null
    }

    // -----------------------------------------------------------------------
    // D3: inference matches on the HOST, exactly or as a dot-suffix. A `contains('openai')` rule
    // would ship a credential to https://evil.example.com/openai/v1.
    // -----------------------------------------------------------------------

    def 'inferProviderFromUrl matches a well-known host and nothing else'() {
        expect:
        AgentConfig.inferProviderFromUrl(ENDPOINT) == PROVIDER

        where:
        ENDPOINT                                         || PROVIDER
        'https://api.openai.com/v1'                      || 'openai'
        'https://api.anthropic.com'                      || 'anthropic'
        'https://openrouter.ai/api/v1'                   || 'openrouter'
        'https://api.mistral.ai/v1'                      || 'mistral'
        and: 'a dot-suffix is the same provider (a regional or versioned subdomain)'
        'https://eu.openrouter.ai/api/v1'                || 'openrouter'
        'https://a.b.api.openai.com/v1'                  || 'openai'
        and: 'case, port, path, query and userinfo are all ignored'
        'https://API.OpenAI.COM/v1'                      || 'openai'
        'https://api.openai.com:8443/v1'                 || 'openai'
        'https://user:pw@api.openai.com/v1?x=1'          || 'openai'
        '  https://api.openai.com/v1  '                  || 'openai'
        and: 'the absolute-FQDN spelling cannot slip past an exact match'
        'https://api.openai.com./v1'                     || 'openai'
        and: '-- REJECTED -- the provider name in the PATH is not the host'
        'https://evil.example.com/openai/v1'             || null
        'https://gateway.corp/openrouter.ai/v1'          || null
        and: '-- REJECTED -- a SUBSTRING of the host is not the host'
        'https://notopenrouter.ai/api/v1'                || null
        'https://openrouter.ai.evil.example/api/v1'      || null
        'https://api.openai.com.evil.example/v1'         || null
        and: '-- REJECTED -- a lookalike that only shares a suffix boundary the table does not own'
        'https://openai.com/v1'                          || null
        'https://myapi.mistral.ai.co/v1'                 || null
        and: '-- REJECTED -- the host in the USERINFO is not the host'
        'https://api.openai.com@evil.example/v1'         || null
        and: 'no host, no answer -- never a guess'
        'http://localhost:8000/v1'                       || null
        'not a url at all'                               || null
        '/v1/chat'                                       || null
        ''                                               || null
        null                                             || null
    }

    // -----------------------------------------------------------------------
    // D2 tier 3, and the gate on it: `<PROVIDER>_API_KEY` is a variable exported for a PROVIDER,
    // not for this endpoint, and a runner installs what it is handed ahead of anything it could
    // resolve itself -- so it travels only when the endpoint agrees.
    // -----------------------------------------------------------------------

    def 'apiKeyFor resolves the provider tier only for an endpoint that provider owns'() {
        expect:
        new AgentConfig(OPTS, ENV).apiKeyFor(MODEL) == KEY

        where:
        ENV                                                                       | OPTS                                                            | MODEL                          || KEY
        // -- the neutral tiers are never gated
        [:]                                                                       | [apiKey: 'sk-cfg', baseUrl: 'https://api.anthropic.com/v1']     | 'openai/gpt-5-mini'            || 'sk-cfg'
        [NXF_AGENT_API_KEY: 'sk-nxf']                                             | [baseUrl: 'https://api.anthropic.com/v1']                       | 'openai/gpt-5-mini'            || 'sk-nxf'
        // -- tier 3, no endpoint at all: the request goes to the provider's own default
        [OPENAI_API_KEY: 'sk-oai']                                                | [:]                                                             | 'openai/gpt-5-mini'            || 'sk-oai'
        [ANTHROPIC_API_KEY: 'sk-ant']                                             | [:]                                                             | 'anthropic/claude-sonnet-4'    || 'sk-ant'
        [MISTRAL_API_KEY: 'sk-mis']                                               | [:]                                                             | 'mistral/mistral-large'        || 'sk-mis'
        [OPENROUTER_API_KEY: 'sk-or']                                             | [:]                                                             | 'openrouter/openai/gpt-4o'     || 'sk-or'
        [AZURE_OPENAI_API_KEY: 'sk-az']                                           | [:]                                                             | 'azure/gpt-4o'                 || 'sk-az'
        [GOOGLE_API_KEY: 'sk-goo']                                                | [:]                                                             | 'gemini/gemini-2.0-flash'      || 'sk-goo'
        [GEMINI_API_KEY: 'sk-gem']                                                | [:]                                                             | 'google/gemini-2.0-flash'      || 'sk-gem'
        // -- tier 3 to the provider's OWN host
        [ANTHROPIC_API_KEY: 'sk-ant']                                             | [baseUrl: 'https://api.anthropic.com/v1']                       | 'anthropic/claude-sonnet-4'    || 'sk-ant'
        and: 'the headline D1 case: openai PROTOCOL, openrouter CREDENTIAL, inferred from the host'
        [OPENROUTER_API_KEY: 'sk-or', OPENAI_API_KEY: 'sk-oai']                   | [baseUrl: 'https://openrouter.ai/api/v1']                       | 'openai/gpt-4o'                || 'sk-or'
        and: 'tier 3 to the endpoint that same namespace supplied'
        [OPENAI_API_KEY: 'sk-oai', OPENAI_BASE_URL: 'http://mirror:8000/v1']      | [:]                                                             | 'openai/gpt-5-mini'            || 'sk-oai'
        and: 'tier 3 to an unrecognized gateway, VOUCHED FOR by an explicit apiProvider'
        [OPENAI_API_KEY: 'sk-oai']                                                | [baseUrl: 'https://gw.corp/v1', apiProvider: 'openai']          | 'openai/gpt-5-mini'            || 'sk-oai'
        and: '-- WITHHELD -- an unrecognized gateway nobody vouched for (the pre-existing misroute)'
        [OPENAI_API_KEY: 'sk-oai']                                                | [baseUrl: 'https://gw.corp/v1']                                 | 'openai/gpt-5-mini'            || null
        and: '-- WITHHELD -- another provider\'s own host, even with an explicit apiProvider'
        [OPENAI_API_KEY: 'sk-oai']                                                | [baseUrl: 'https://api.anthropic.com/v1', apiProvider: 'openai'] | 'openai/gpt-5-mini'           || null
        and: '-- WITHHELD -- OPENAI_BASE_URL pointed at somebody else closes the same hole'
        [OPENAI_API_KEY: 'sk-oai', OPENAI_BASE_URL: 'https://openrouter.ai/api/v1'] | [:]                                                           | 'openai/gpt-5-mini'            || null
        and: '-- NOTHING TO READ -- the wrong provider\'s variable is never consulted'
        [OPENAI_API_KEY: 'sk-oai']                                                | [:]                                                             | 'anthropic/claude-sonnet-4'    || null
        [ANTHROPIC_API_KEY: 'sk-ant']                                             | [:]                                                             | 'openai/gpt-5-mini'            || null
        and: '-- NOTHING TO READ -- an unknown namespace names no variable at all'
        [STRIPE_API_KEY: 'sk-stripe']                                             | [:]                                                             | 'stripe/whatever'              || null
        and: '-- NOTHING TO READ -- no prefix means no provider'
        [OPENAI_API_KEY: 'sk-oai']                                                | [:]                                                             | 'gpt-5-mini'                   || null
        [OPENAI_API_KEY: 'sk-oai']                                                | [:]                                                             | null                           || null
        [:]                                                                       | [:]                                                             | 'openai/gpt-5-mini'            || null
    }

    def 'baseUrlFor resolves the provider tier ungated -- an endpoint is not a secret'() {
        expect:
        new AgentConfig(OPTS, ENV).baseUrlFor(MODEL) == ENDPOINT

        where:
        ENV                                                                  | OPTS                       | MODEL                       || ENDPOINT
        [:]                                                                  | [baseUrl: 'http://cfg/v1'] | 'openai/gpt-4o'             || 'http://cfg/v1'
        [NXF_AGENT_BASE_URL: 'http://nxf/v1']                                | [:]                        | 'anthropic/claude-sonnet-4' || 'http://nxf/v1'
        [NXF_AGENT_BASE_URL: 'http://nxf/v1', OPENAI_BASE_URL: 'http://o/v1']| [:]                        | 'openai/gpt-4o'             || 'http://nxf/v1'
        and: 'tier 3, per provider, on the three variables that actually exist'
        [OPENAI_BASE_URL: 'http://oai/v1']                                   | [:]                        | 'openai/gpt-4o'             || 'http://oai/v1'
        [ANTHROPIC_BASE_URL: 'http://ant/v1']                                | [:]                        | 'anthropic/claude-sonnet-4' || 'http://ant/v1'
        [AZURE_OPENAI_ENDPOINT: 'https://x.openai.azure.com']                | [:]                        | 'azure/gpt-4o'              || 'https://x.openai.azure.com'
        and: 'the variable is scoped to ITS provider, so another model does not pick it up'
        [OPENAI_BASE_URL: 'http://oai/v1']                                   | [:]                        | 'anthropic/claude-sonnet-4' || null
        and: 'an explicit apiProvider redirects which variable is read'
        [OPENAI_BASE_URL: 'http://oai/v1']                                   | [apiProvider: 'anthropic'] | 'openai/gpt-4o'             || null
        [ANTHROPIC_BASE_URL: 'http://ant/v1']                                | [apiProvider: 'anthropic'] | 'openai/gpt-4o'             || 'http://ant/v1'
        and: 'no invented conventions: these variables are read by nobody'
        [MISTRAL_BASE_URL: 'http://mis/v1']                                  | [:]                        | 'mistral/mistral-large'     || null
        [OPENROUTER_BASE_URL: 'http://or/v1']                                | [:]                        | 'openrouter/x'              || null
        [GOOGLE_BASE_URL: 'http://goo/v1']                                   | [:]                        | 'google/gemini-2.0-flash'   || null
        and: 'nothing set means "use the provider default"'
        [:]                                                                  | [:]                        | 'openai/gpt-4o'             || null
        [OPENAI_BASE_URL: 'http://oai/v1']                                   | [:]                        | null                        || null
    }

    def 'a withheld provider credential is named, with both remedies'() {
        given: 'silence would leave an opaque 401 as the only evidence that a key was resolved and dropped'
        def config = new AgentConfig([baseUrl: 'https://gw.corp/v1'], [OPENAI_API_KEY: 'sk-oai'])

        when:
        def resolved = config.apiKeyFor('openai/gpt-5-mini')

        then:
        resolved == null
        and: 'the WARN names the variable that was read, the endpoint, and the two ways to proceed'
        capture.toString().contains('Not using the OPENAI_API_KEY credential')
        capture.toString().contains('https://gw.corp/v1')
        capture.toString().contains('`agent.apiKey`')
        capture.toString().contains("`agent.apiProvider = 'openai'`")
    }

    def 'an explicit apiProvider contradicting a well-known endpoint host is called out at build time'() {
        when: 'the endpoint is Anthropic\'s own API but the config claims an openai credential namespace'
        // far likelier a mistake than an intention -- and the mistake would ship a credential to a
        // third party, so it is said ONCE per agent here rather than only when a key is withheld
        new AgentConfig([apiProvider: 'openai', baseUrl: 'https://api.anthropic.com/v1'], [:])

        then:
        capture.toString().contains("agent.apiProvider = 'openai'")
        capture.toString().contains('is a known `anthropic` endpoint')

        when: 'the two agree, or the host is not one the table recognizes'
        def before = capture.toString().length()
        new AgentConfig([apiProvider: 'anthropic', baseUrl: 'https://api.anthropic.com/v1'], [:])
        new AgentConfig([apiProvider: 'openai', baseUrl: 'https://gw.corp/v1'], [:])
        new AgentConfig([apiProvider: 'openai'], [:])

        then: 'there is nothing to contradict, so nothing is said'
        !capture.toString().substring(before).contains('is a known')
    }

    def 'missingCredentialHint names the variables the ladder ACTUALLY consults'() {
        expect: 'a message naming OPENAI_API_KEY for an anthropic run is worse than no message'
        AgentConfig.missingCredentialHint('anthropic').contains('`NXF_AGENT_API_KEY` or `ANTHROPIC_API_KEY`')
        AgentConfig.missingCredentialHint('anthropic').contains('`agent.apiKey`')

        and: 'a provider with an alias pair names both, in resolution order'
        AgentConfig.missingCredentialHint('gemini').contains('`NXF_AGENT_API_KEY` or `GEMINI_API_KEY` or `GOOGLE_API_KEY`')

        and: 'a known provider needs no `agent.apiProvider` advice -- the namespace is already right'
        !AgentConfig.missingCredentialHint('openai').contains('agent.apiProvider')

        and: 'while an unresolvable namespace has only the neutral variable, and IS told to name one'
        AgentConfig.missingCredentialHint('stripe').contains('`NXF_AGENT_API_KEY` environment variable')
        AgentConfig.missingCredentialHint('stripe').contains('`agent.apiProvider`')
        AgentConfig.missingCredentialHint(null).contains('`agent.apiProvider`')
    }

    def 'with no endpoint the provider tier travels only when it IS the model prefix provider'() {
        given: 'no endpoint does not mean "the resolved provider\'s default". The runner dials the'
        // default endpoint of the MODEL PREFIX provider -- the langchain4j client hardcodes
        // https://api.openai.com/v1 -- so `agent.apiProvider = 'openrouter'` with `openai/…` and no
        // `agent.baseUrl` would ship OPENROUTER_API_KEY to OpenAI. It is withheld instead.
        def misrouted = new AgentConfig([apiProvider: 'openrouter'], [OPENROUTER_API_KEY: 'sk-or'])

        expect:
        misrouted.apiKeyFor('openai/gpt-4o') == null
        misrouted.baseUrlFor('openai/gpt-4o') == null
        and: 'and it is a WITHHELD credential, not an absent one'
        misrouted.credentialWithheldFor('openai/gpt-4o')

        and: 'the warn names the variable, why it is not sent, and the endpoint option that fixes it'
        capture.toString().contains('Not using the OPENROUTER_API_KEY credential')
        capture.toString().contains('no endpoint resolved')
        capture.toString().contains('`agent.baseUrl`')

        when: 'the provider IS the prefix, so the default endpoint dialled is that provider\'s own'
        def aligned = new AgentConfig([apiProvider: 'openrouter'], [OPENROUTER_API_KEY: 'sk-or'])

        then:
        aligned.apiKeyFor('openrouter/openai/gpt-4o') == 'sk-or'
        !aligned.credentialWithheldFor('openrouter/openai/gpt-4o')

        when: 'naming the endpoint is what unblocks the redirected namespace'
        def routed = new AgentConfig([apiProvider: 'openrouter', baseUrl: 'https://openrouter.ai/api/v1'],
            [OPENROUTER_API_KEY: 'sk-or'])

        then:
        routed.apiKeyFor('openai/gpt-4o') == 'sk-or'
        !routed.credentialWithheldFor('openai/gpt-4o')
    }

    def 'credentialWithheldFor separates "resolved but refused" from "nothing resolved"'() {
        given: 'the two are indistinguishable through apiKeyFor, which answers null to both -- yet'
        // only the first is a misconfiguration a runner can name, and neither may become the D8
        // placeholder. See AgentRunnerRequest.credential().
        expect:
        new AgentConfig(OPTS, ENV).apiKeyFor(MODEL) == null
        new AgentConfig(OPTS, ENV).credentialWithheldFor(MODEL) == WITHHELD

        where:
        ENV                                                                        | OPTS                                                             | MODEL               || WITHHELD
        // -- RESOLVED AND REFUSED: the key is right there, and the endpoint is not its owner's
        [OPENAI_API_KEY: 'sk-oai']                                                 | [baseUrl: 'https://gw.corp/v1']                                  | 'openai/gpt-5-mini' || true
        [OPENAI_API_KEY: 'sk-oai']                                                 | [baseUrl: 'https://api.anthropic.com/v1', apiProvider: 'openai'] | 'openai/gpt-5-mini' || true
        [OPENAI_API_KEY: 'sk-oai', OPENAI_BASE_URL: 'https://openrouter.ai/api/v1']| [:]                                                              | 'openai/gpt-5-mini' || true
        [OPENROUTER_API_KEY: 'sk-or']                                              | [apiProvider: 'openrouter']                                      | 'openai/gpt-4o'     || true
        and: 'a model id with no provider prefix names no default endpoint either, so nothing is sent'
        [OPENAI_API_KEY: 'sk-oai']                                                 | [apiProvider: 'openai']                                          | 'gpt-5-mini'        || true
        and: '-- NOTHING RESOLVED: no variable was set for the namespace at all'
        [:]                                                                        | [baseUrl: 'https://gw.corp/v1']                                  | 'openai/gpt-5-mini' || false
        [ANTHROPIC_API_KEY: 'sk-ant']                                              | [baseUrl: 'https://gw.corp/v1']                                  | 'openai/gpt-5-mini' || false
        [STRIPE_API_KEY: 'sk-stripe']                                              | [:]                                                              | 'stripe/whatever'   || false
    }

    def 'a resolved or neutral credential is never reported as withheld'() {
        expect: 'tiers 1 and 2 are not gated at all, so there is nothing to withhold'
        !new AgentConfig([apiKey: 'sk-cfg', baseUrl: 'https://gw.corp/v1'], [OPENAI_API_KEY: 'sk-oai'])
            .credentialWithheldFor('openai/gpt-5-mini')
        !new AgentConfig([baseUrl: 'https://gw.corp/v1'], [NXF_AGENT_API_KEY: 'sk-nxf', OPENAI_API_KEY: 'sk-oai'])
            .credentialWithheldFor('openai/gpt-5-mini')

        and: 'nor is a provider tier the gate lets through'
        !new AgentConfig([:], [OPENAI_API_KEY: 'sk-oai']).credentialWithheldFor('openai/gpt-5-mini')
        !new AgentConfig([baseUrl: 'https://api.openai.com/v1'], [OPENAI_API_KEY: 'sk-oai'])
            .credentialWithheldFor('openai/gpt-5-mini')
    }

    def 'the D3 inference that outranks the model prefix is logged, naming both'() {
        given: 'inference is invisible magic otherwise: which variable was read cannot be told from'
        // the config, and reading the WRONG one is a credential misroute. Logged at debug on the
        // resolution path, so it appears exactly when a credential was actually resolved.
        when: 'the endpoint host names the provider and the prefix names a different one'
        def fromHost = new AgentConfig([baseUrl: 'https://openrouter.ai/api/v1'], [OPENROUTER_API_KEY: 'sk-or'])
            .apiKeyFor('openai/gpt-4o')

        then:
        fromHost == 'sk-or'
        capture.toString().contains('Resolved API provider `openrouter` for agent model `openai/gpt-4o`')
        capture.toString().contains('from the endpoint https://openrouter.ai/api/v1')
        capture.toString().contains('(the model prefix is `openai`)')

        when: 'an explicit option is what redirected the namespace, the line says so instead'
        def before = capture.toString().length()
        def fromOption = new AgentConfig([apiProvider: 'openrouter', baseUrl: 'https://openrouter.ai/api/v1'],
            [OPENROUTER_API_KEY: 'sk-or']).apiKeyFor('openai/gpt-4o')

        then:
        fromOption == 'sk-or'
        capture.toString().substring(before).contains('from `agent.apiProvider`')

        when: 'the provider and the prefix agree, there is no surprise to report'
        def quiet = capture.toString().length()
        new AgentConfig([:], [OPENAI_API_KEY: 'sk-oai']).apiKeyFor('openai/gpt-5-mini')

        then:
        !capture.toString().substring(quiet).contains('Resolved API provider')
    }

    def 'a typo inside the nested rpc scope is reported'() {
        given: 'the whole point of declaring `rpc` as a nested ConfigScope rather than a Map option:'
        // a Map option short-circuits ConfigValidator.isMapOption and the sub-map is never walked,
        // so every key under it -- typo or not -- was silently accepted
        when:
        new ConfigValidator().validate([agent: [rpc: [remoteHostt: 'driver.internal']]])
        then:
        capture.toString().contains("Unrecognized config option 'agent.rpc.remoteHostt'")
    }


    // NOTE: the positive diagnostic for an agent BODY directive name written in the config
    // (`agent.model`, `agent.maxIterations`) lives in ConfigValidatorTest. `log.warn1` caches
    // by message for the whole JVM, so asserting the same warning twice across specs is
    // inherently order-dependent -- see the existing `wokDir`/`wokDir2` split there.
}
