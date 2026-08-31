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

import nextflow.Global
import nextflow.Session
import nextflow.SysEnv
import nextflow.exception.ConfigParseException
import nextflow.exception.IllegalConfigException
import nextflow.exception.ScriptRuntimeException
import nextflow.processor.ConfigList
import nextflow.script.AgentBuilder.AgentInput
import nextflow.script.AgentBuilder.AgentOutput
import nextflow.script.AgentDef
import nextflow.script.BaseScript
import nextflow.script.PromptDef
import nextflow.script.ScriptBinding
import org.junit.Rule
import spock.lang.Timeout
import test.Dsl2Spec
import test.MockSession
import test.OutputCapture

/**
 * The `agent` config scope resolves task directives with the SAME semantics as the
 * `process` scope (selectors, precedence, `ext` merge, repeatable directives) while
 * staying fully independent from it.
 *
 * The selector matchers themselves ({@code matchesSelector}/{@code matchesLabels}) and
 * the repeat/ext machinery are reused verbatim from
 * {@link nextflow.script.dsl.ProcessConfigBuilder} and are pinned by
 * {@code ProcessConfigBuilderTest}; this spec pins the AGENT wiring of that machinery.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Timeout(30)
class AgentConfigSelectorTest extends Dsl2Spec {

    @Rule
    OutputCapture capture = new OutputCapture()

    def cleanup() {
        AgentRunnerProvider.testRunner = null
    }

    /**
     * Spin a MockSession (MockExecutorFactory) and make it the global session, adding the
     * configuration that CONTAINERIZES a canonical agent task -- an enabled container engine plus an
     * `agent.container` image, which every launch-spec runner now requires because its launch command
     * is built from paths that exist only inside the runner image. This spec is about selector
     * RESOLUTION, so the containerization is boilerplate here; the tests that deliberately exercise a
     * missing or disabled image use {@link #bareSession} instead.
     *
     * <p>An `agent.container` already present in the PLAIN scope is preserved, and a selector-provided
     * one still wins over the injected plain value by the normal precedence ladder.
     */
    private Session newSession(Map config = null) {
        final Map merged = new LinkedHashMap(config ?: [:])
        merged.docker = [enabled: true]
        final Map agentScope = new LinkedHashMap((Map) (merged.agent ?: [:]))
        agentScope.container = agentScope.container ?: 'agent-image:test'
        merged.agent = agentScope
        return bareSession(merged)
    }

    /** Spin a MockSession with EXACTLY the given config (no containerization defaults). */
    private Session bareSession(Map config = null) {
        def session = config ? new MockSession(config) : new MockSession()
        session.setBinding(new ScriptBinding())
        session.init(null, null, null, null)
        session.start()
        Global.session = session
        return session
    }

    private AgentDef newAgent(Map directives = [model: 'openai/gpt-4o']) {
        final owner = Mock(BaseScript) { getBinding() >> new ScriptBinding() }
        return new AgentDef(owner, 'qa', directives as Map<String,Object>,
            [new AgentInput('q', String)], [new AgentOutput('answer', String)],
            new PromptDef({ -> 'Q' }, 'Q'))
    }

    /**
     * A launch-spec runner. With no argument it declares NO image of its own, which is the shape
     * every containerization spec below assumes; pass a coordinate for the pi-shaped runner that
     * generates one from its own VERSION.
     */
    private static AgentRunner canonicalRunner(String defaultContainer = null) {
        return new AgentRunner() {
            @Override String getName() { 'canonical-test' }
            @Override AgentLaunchSpec getLaunchSpec() {
                new AgentLaunchSpec(['/agent-rpc'], ['node'])
            }
            @Override String getDefaultContainer() { defaultContainer }
            @Override String run(AgentRunnerRequest req) { 'x' }
        }
    }

    def 'a withName selector beats the plain agent scope and only matches its target'() {
        given:
        AgentRunnerProvider.testRunner = canonicalRunner()

        when:
        newSession([agent: [cpus: 1, 'withName:qa': [cpus: 2]]])
        then:
        newAgent().buildAgentTask(['hello']).config.cpus == 2

        when: 'the selector names a different agent'
        newSession([agent: [cpus: 1, 'withName:other': [cpus: 2]]])
        then:
        newAgent().buildAgentTask(['hello']).config.cpus == 1
    }

    def 'a withLabel selector matches an agent declaring that label'() {
        given:
        AgentRunnerProvider.testRunner = canonicalRunner()

        when: 'the agent declares the label'
        newSession([agent: [cpus: 1, 'withLabel:reasoning': [cpus: 4]]])
        then:
        newAgent([model: 'openai/gpt-4o', label: ['reasoning']]).buildAgentTask(['hello']).config.cpus == 4

        when: 'the agent declares no label'
        newSession([agent: [cpus: 1, 'withLabel:reasoning': [cpus: 4]]])
        then:
        newAgent().buildAgentTask(['hello']).config.cpus == 1

        when: 'a negated label selector faces a label-less agent'
        newSession([agent: [cpus: 1, 'withLabel:!reasoning': [cpus: 4]]])
        then:
        newAgent().buildAgentTask(['hello']).config.cpus == 4
    }

    def 'the selector precedence ladder resolves fully-qualified over alias, base, label and plain'() {
        given:
        AgentRunnerProvider.testRunner = canonicalRunner()
        def scope = [
            cpus: 1,
            'withLabel:reasoning': [cpus: 2],
            'withName:qa': [cpus: 3],
            'withName:WF:qa': [cpus: 4],
        ]
        def agent = { -> (AgentDef) newAgent([model: 'openai/gpt-4o', label: ['reasoning']]).cloneWithName('WF:qa') }

        when:
        newSession([agent: scope])
        then: 'the fully-qualified rule wins'
        agent().buildAgentTask(['hello']).config.cpus == 4

        when: 'the fully-qualified rule is removed'
        newSession([agent: scope.findAll { it.key != 'withName:WF:qa' }])
        then: 'the base-name rule wins'
        agent().buildAgentTask(['hello']).config.cpus == 3

        when: 'the name rules are removed'
        newSession([agent: [cpus: 1, 'withLabel:reasoning': [cpus: 2]]])
        then: 'the label rule wins'
        agent().buildAgentTask(['hello']).config.cpus == 2

        when: 'no selector is left'
        newSession([agent: [cpus: 1]])
        then: 'the plain scope applies'
        agent().buildAgentTask(['hello']).config.cpus == 1
    }

    def 'an aliased include still matches the declared name and the alias'() {
        given:
        AgentRunnerProvider.testRunner = canonicalRunner()

        when: 'the selector names the DECLARED agent (baseName is preserved by cloneWithName)'
        newSession([agent: ['withName:qa': [cpus: 3]]])
        then:
        newAgent().cloneWithName('reviewer').buildAgentTask(['hello']).config.cpus == 3

        when: 'the selector names the alias'
        newSession([agent: ['withName:reviewer': [cpus: 5]]])
        then:
        newAgent().cloneWithName('reviewer').buildAgentTask(['hello']).config.cpus == 5
    }

    def 'ext maps merge instead of replacing'() {
        given:
        AgentRunnerProvider.testRunner = canonicalRunner()

        when: 'plain ext and a selector ext declare different keys'
        newSession([agent: [ext: [args: '--verbose'], 'withName:qa': [ext: [opts: '--fast']]]])
        def config = newAgent().buildAgentTask(['hello']).config
        then: 'both survive'
        config.ext == [args: '--verbose', opts: '--fast']

        when: 'the selector overwrites the same key'
        newSession([agent: [ext: [args: '--verbose'], 'withName:qa': [ext: [args: '--fast']]]])
        then: 'the selector wins and the plain value does not come back'
        newAgent().buildAgentTask(['hello']).config.ext == [args: '--fast']
    }

    def 'the label directive goes through the repeatable-directive path'() {
        given:
        AgentRunnerProvider.testRunner = canonicalRunner()

        when: 'only a selector sets the label'
        newSession([agent: ['withName:qa': [label: 'x']]])
        def config = newAgent().buildAgentTask(['hello']).config
        then: 'the value is a validated ConfigList, as for a process'
        config.label instanceof ConfigList
        config.getLabels() == ['x']

        when: 'the body declares a label and a matching withLabel selector sets another'
        newSession([agent: ['withLabel:reasoning': [label: 'fast']]])
        def config2 = newAgent([model: 'openai/gpt-4o', label: ['reasoning']]).buildAgentTask(['hello']).config
        then: 'exact `process` parity: ProcessConfigBuilder.putWithRepeat REPLACES the declared list,'
        // and the selector still matched, because applyConfig reads the DECLARED labels on entry
        config2.getLabels() == ['fast']

        when: 'a body label is declared and no selector touches `label`'
        newSession([agent: [cpus: 1]])
        then: 'the declared labels reach the config'
        newAgent([model: 'openai/gpt-4o', label: ['reasoning', 'fast']]).buildAgentTask(['hello']).config.getLabels() == ['reasoning', 'fast']

        when: 'the declared label is not a valid label'
        newSession([agent: [cpus: 1]])
        newAgent([model: 'openai/gpt-4o', label: ['a-b']]).buildAgentTask(['hello'])
        then: 'the diagnostic names the AGENT: `label` is the one directive method agents share'
        def error = thrown(IllegalConfigException)
        error.message.startsWith('Not a valid agent label: a-b')
    }

    def 'agent-only options never reach the task config and never warn'() {
        given:
        AgentRunnerProvider.testRunner = canonicalRunner()
        // `rpc` is a NESTED ConfigScope field, not a @ConfigOption: it must still be excluded from
        // the directive axis, or it lands in every agent's config as a phantom `rpc` directive
        // (silently in the plain scope, as `Unknown directive` inside a selector)
        def agentOnly = ['runner', 'model', 'apiProvider', 'apiKey', 'baseUrl', 'maxIterations', 'requestTimeout',
                         'maxToolOutputInlineSize', 'trace', 'rpc']

        when: 'the agent-only options are set in the plain scope'
        newSession([agent: [runner: 'test', model: 'm', apiProvider: 'openai', apiKey: 'sk-x', baseUrl: 'http://x/v1', maxIterations: 7, rpc: [port: 1234], cpus: 1]])
        def config = newAgent().buildAgentTask(['hello']).config
        then:
        agentOnly.every { !config.containsKey(it) }
        config.cpus == 1
        !capture.toString().contains('Unknown directive')

        when: 'the same options are set inside a selector'
        newSession([agent: ['withName:qa': [runner: 'test', model: 'm', apiProvider: 'openai', apiKey: 'sk-x', baseUrl: 'http://x/v1', maxIterations: 7, rpc: [port: 1234], cpus: 2]]])
        def config2 = newAgent().buildAgentTask(['hello']).config
        then:
        agentOnly.every { !config2.containsKey(it) }
        config2.cpus == 2
        !capture.toString().contains('Unknown directive')
    }

    def 'an unknown directive is reported with the agent noun'() {
        given:
        AgentRunnerProvider.testRunner = canonicalRunner()
        newSession([agent: ['withName:qa': [fooBar: 1]]])

        when:
        newAgent().buildAgentTask(['hello'])

        then:
        capture.toString().contains('Unknown directive `fooBar` for agent `qa`')
    }

    def 'the local executor default does not shadow a selector'() {
        given:
        AgentRunnerProvider.testRunner = canonicalRunner()

        when: 'no executor is configured'
        newSession([agent: [cpus: 1]])
        then:
        newAgent().buildAgentTask(['hello']).config.executor == AgentConfig.DEFAULT_EXECUTOR

        // an offloaded agent must also declare a container, and -- because its container is
        // launched off the driver host -- an address the container can reach the driver on
        when: 'a selector sets the executor'
        newSession([agent: [
            rpc: [remoteHost: 'driver.internal'],
            'withName:qa': [executor: 'k8s', container: 'agent-image:1'] ]])
        then:
        newAgent().buildAgentTask(['hello']).config.executor == 'k8s'
    }

    def 'an explicitly declared agent.container resolves through the ladder'() {
        given:
        AgentRunnerProvider.testRunner = canonicalRunner()

        when: 'the plain scope sets a container'
        bareSession([docker: [enabled: true], agent: [container: 'x']])
        then:
        newAgent().buildAgentTask(['hello']).config.container == 'x'

        when: 'a selector sets a container'
        bareSession([docker: [enabled: true], agent: ['withName:qa': [container: 'y']]])
        then:
        newAgent().buildAgentTask(['hello']).config.container == 'y'

        when: 'a selector overrides the plain scope'
        bareSession([docker: [enabled: true], agent: [container: 'x', 'withName:qa': [container: 'y']]])
        then:
        newAgent().buildAgentTask(['hello']).config.container == 'y'
    }

    // -- a runner whose runtime lives IN an image knows which image that is (nf-agent-pi generates
    //    the coordinate from its own VERSION), so `agent.container` is optional for it. The value
    //    lands in the config exactly where an explicit one would, so nothing downstream -- the
    //    containerization guard, the per-task re-check, the container fingerprint in the task hash
    //    -- has to know it was defaulted.
    def 'a runner that declares an image of its own defaults agent.container to it'() {
        given:
        AgentRunnerProvider.testRunner = canonicalRunner('registry.io/nf-agent-pi:1.2.3')

        when: 'nothing in the config declares a container'
        bareSession([docker: [enabled: true], agent: [:]])
        then:
        newAgent().buildAgentTask(['hello']).config.container == 'registry.io/nf-agent-pi:1.2.3'

        when: 'and with the agent offloaded to a remote executor'
        bareSession([docker: [enabled: true], agent: [executor: 'k8s', rpc: [remoteHost: 'driver.internal']]])
        then:
        newAgent().buildAgentTask(['hello']).config.container == 'registry.io/nf-agent-pi:1.2.3'
    }

    def 'an explicit agent.container beats the runner image everywhere on the ladder'() {
        given:
        AgentRunnerProvider.testRunner = canonicalRunner('registry.io/nf-agent-pi:1.2.3')

        when: 'the plain scope declares one'
        bareSession([docker: [enabled: true], agent: [container: 'mine:1']])
        then:
        newAgent().buildAgentTask(['hello']).config.container == 'mine:1'

        when: 'a `withLabel:` selector declares one -- the default must not pre-empt the ladder'
        bareSession([docker: [enabled: true], agent: ['withLabel:reasoning': [container: 'mine:2']]])
        then:
        newAgent([model: 'openai/gpt-4o', label: ['reasoning']]).buildAgentTask(['hello']).config.container == 'mine:2'
    }

    def 'the explicit container opt-out survives a runner that declares an image'() {
        given: '`agent.container = false` is present-with-value-false, not absent'
        AgentRunnerProvider.testRunner = canonicalRunner('registry.io/nf-agent-pi:1.2.3')
        bareSession([docker: [enabled: true], agent: [container: false]])

        when:
        newAgent().buildAgentTask(['hello'])

        then: 'the opt-out keeps meaning "no container", so the ORIGINAL message still fires'
        def e = thrown(ScriptRuntimeException)
        e.message.contains('must declare a container')
    }

    // -- the most common first run of a pi agent: the plugin is installed, no `agent.container` is
    //    declared and no container engine is enabled. With the image now defaulted, the `!hasContainer`
    //    branch is unreachable and this is the message the user gets, so it must not tell them they
    //    declared something they did not.
    def 'a defaulted image with no engine enabled is rejected without blaming the user for it'() {
        given:
        AgentRunnerProvider.testRunner = canonicalRunner('registry.io/nf-agent-pi:1.2.3')
        bareSession([agent: [:]])

        when:
        newAgent().buildAgentTask(['hello'])

        then:
        def e = thrown(ScriptRuntimeException)
        e.message.contains('would not run it in a container')
        e.message.contains('docker.enabled')
        and: 'the user declared no `agent.container`, so the message must not say they did'
        !e.message.contains('declares `agent.container`')
    }

    def 'a canonical agent without a container fails fast on EVERY executor'() {
        given: 'the launch command is built from paths that exist only inside the runner image'
        AgentRunnerProvider.testRunner = canonicalRunner()

        when: 'nothing is configured -- and this runner declares no image of its own'
        bareSession([docker: [enabled: true], agent: [:]])
        newAgent().buildAgentTask(['hello'])

        then:
        def local = thrown(ScriptRuntimeException)
        local.message.contains('must declare a container')

        when: 'the agent is offloaded to a remote executor with no image'
        bareSession([docker: [enabled: true], agent: [executor: 'k8s']])
        newAgent().buildAgentTask(['hello'])

        then:
        def offloaded = thrown(ScriptRuntimeException)
        offloaded.message.contains('must declare a container')

        when: 'the container is explicitly disabled (present-with-value-false)'
        bareSession([docker: [enabled: true], agent: [executor: 'k8s', container: false]])
        newAgent().buildAgentTask(['hello'])

        then: 'the opt-out is rejected exactly like an absent image'
        def disabledError = thrown(ScriptRuntimeException)
        disabledError.message.contains('must declare a container')

        // `rpc.remoteHost` joins the image as a requirement once the executor is remote: the
        // container is launched off the driver host, so no engine alias can name the driver
        when: 'a container is declared'
        bareSession([docker: [enabled: true], agent: [
            executor: 'k8s', container: 'agent-image:1', rpc: [remoteHost: 'driver.internal'] ]])
        then:
        newAgent().buildAgentTask(['hello']).config.container == 'agent-image:1'
    }

    def 'a canonical agent with an image but no container engine fails fast'() {
        given: 'a non-container-native executor with every engine disabled would run the in-image paths on the host'
        AgentRunnerProvider.testRunner = canonicalRunner()
        bareSession([agent: [container: 'agent-image:1']])

        when:
        newAgent().buildAgentTask(['hello'])

        then:
        def error = thrown(ScriptRuntimeException)
        error.message.contains('would not run it in a container')
        error.message.contains('docker.enabled')
    }

    def 'a process selector never configures an agent'() {
        given:
        AgentRunnerProvider.testRunner = canonicalRunner()
        newSession([process: ['withName:qa': [cpus: 9]], agent: [cpus: 1]])

        expect:
        newAgent().buildAgentTask(['hello']).config.cpus == 1
    }

    def 'a plain value beats a matching selector when the selector value equals the built-in default'() {
        given: 'the inherited `process` quirk: applyConfigDefaults re-applies a plain value when the'
        // current value still EQUALS ProcessConfig.DEFAULT_CONFIG (ProcessConfigBuilder.applyConfigDefaults),
        // and DEFAULT_CONFIG.maxRetries == 1. This is exact `process` parity, deliberately not fixed.
        AgentRunnerProvider.testRunner = canonicalRunner()
        newSession([agent: [maxRetries: 5, 'withName:qa': [maxRetries: 1]]])

        expect:
        newAgent().buildAgentTask(['hello']).config.maxRetries == 5
    }

    def 'the agent-only axis is selector-resolved end to end'() {
        given:
        AgentRunnerProvider.testRunner = canonicalRunner()

        when:
        newSession([agent: [model: 'm0', 'withName:qa': [model: 'openai/gpt-5', maxIterations: 40]]])
        def config = newAgent([:]).agentConfig()
        then:
        config.model == 'openai/gpt-5'
        config.maxIterations == 40

        when: 'an agent declaring neither `model` nor `maxIterations` is built'
        // the EFFECTIVE resolved values enter the canonical task identity (BodyDef.source),
        // so this reads what buildAgentTask actually consumed
        def defaulted = newAgent([:]).buildAgentTask(['hello']).getTaskBody().source
        then: 'the selector-resolved defaults apply'
        defaulted.contains('agentModel=openai/gpt-5\n')
        defaulted.contains('maxIterations=40\n')

        when: 'the agent declares its own model and maxIterations'
        def declared = newAgent([model: 'openai/gpt-4o', maxIterations: 3]).buildAgentTask(['hello']).getTaskBody().source
        then: 'the body directives still win -- the config options are DEFAULTS on every rung'
        declared.contains('agentModel=openai/gpt-4o\n')
        declared.contains('maxIterations=3\n')
    }

    def 'a selector setting only apiKey leaves the outer baseUrl intact (the per-key merge)'() {
        given: 'this is WHY the endpoint and the credential are FLAT options (design D3):'
        // AgentConfig.copyOptions merges key by key, so a selector overriding one leaf cannot
        // silently drop its sibling. A nested `agent.config { }` container would make the whole
        // container the unit of copy, and a selector setting one leaf would lose the other.
        AgentRunnerProvider.testRunner = canonicalRunner()
        and: 'an empty environment, so an exported OPENAI_* on the dev machine cannot supply a tier'
        SysEnv.push([:])

        when: 'the selector overrides only the credential'
        newSession([agent: [apiKey: 'sk-outer', baseUrl: 'http://outer/v1', 'withName:qa': [apiKey: 'sk-inner']]])
        def merged = newAgent([:]).agentConfig()
        then: 'the selector wins for `apiKey` and the outer `baseUrl` SURVIVES'
        merged.apiKey == 'sk-inner'
        merged.baseUrl == 'http://outer/v1'

        when: 'the selector overrides only the endpoint'
        newSession([agent: [apiKey: 'sk-outer', baseUrl: 'http://outer/v1', 'withName:qa': [baseUrl: 'http://inner/v1']]])
        def swapped = newAgent([:]).agentConfig()
        then: 'symmetrically, the outer credential survives'
        swapped.apiKey == 'sk-outer'
        swapped.baseUrl == 'http://inner/v1'

        when: 'the selector names a DIFFERENT agent'
        newSession([agent: [apiKey: 'sk-outer', baseUrl: 'http://outer/v1', 'withName:other': [apiKey: 'sk-inner']]])
        def untouched = newAgent([:]).agentConfig()
        then: 'the plain scope applies to both keys'
        untouched.apiKey == 'sk-outer'
        untouched.baseUrl == 'http://outer/v1'

        when: 'only a selector declares them (nothing in the plain scope)'
        newSession([agent: ['withName:qa': [apiKey: 'sk-inner', baseUrl: 'http://inner/v1']]])
        def selectorOnly = newAgent([:]).agentConfig()
        then:
        selectorOnly.apiKey == 'sk-inner'
        selectorOnly.baseUrl == 'http://inner/v1'

        cleanup:
        SysEnv.pop()
    }

    def 'a selector setting only apiProvider names the credential namespace for that agent alone'() {
        given: 'the namespace is per-agent like every other agent-only option (design D1): one'
        // pipeline can reach OpenRouter through the openai protocol for one agent and OpenAI proper
        // for another, and `apiProvider` is the only thing that tells the two credentials apart
        AgentRunnerProvider.testRunner = canonicalRunner()
        SysEnv.push([OPENAI_API_KEY: 'sk-openai', OPENROUTER_API_KEY: 'sk-openrouter'])

        when: 'only the selector names the namespace, over a plain-scope endpoint'
        newSession([agent: [baseUrl: 'https://gw.corp/v1', 'withName:qa': [apiProvider: 'openrouter']]])
        def selected = newAgent([:]).agentConfig()
        then: 'the sibling endpoint survives the per-key merge, and the namespace redirects the credential'
        selected.apiProvider == 'openrouter'
        selected.baseUrl == 'https://gw.corp/v1'
        selected.apiKeyFor('openai/gpt-4o') == 'sk-openrouter'

        when: 'the selector names a DIFFERENT agent'
        newSession([agent: [baseUrl: 'https://gw.corp/v1', 'withName:other': [apiProvider: 'openrouter']]])
        def untouched = newAgent([:]).agentConfig()
        then: 'nothing vouches for the gateway, so the ambient provider key is withheld from it'
        untouched.apiProvider == null
        untouched.apiKeyFor('openai/gpt-4o') == null

        cleanup:
        SysEnv.pop()
    }

    def 'a dynamic directive resolves against the task context'() {
        given: 'the nf-lang `isProcessScope` fix makes `task` resolvable in an agent-scope closure;'
        // this pins that the closure survives applyConfig into the agent's ProcessConfigV2
        // and is evaluated per task, not stored as a literal
        AgentRunnerProvider.testRunner = canonicalRunner()
        newSession([agent: ['withName:qa': [cpus: { task.attempt * 2 }]]])

        when:
        def taskConfig = newAgent().buildAgentTask(['hello']).config.createTaskConfig()
        taskConfig.setContext([:])
        then: 'attempt defaults to 1'
        taskConfig.getCpus() == 2
    }

    def 'a malformed selector body is reported as a config parse error'() {
        given:
        AgentRunnerProvider.testRunner = canonicalRunner()
        newSession([agent: ['withName:qa': 'oops']])

        when:
        newAgent().buildAgentTask(['hello'])

        then:
        thrown(ConfigParseException)
    }

    // NOTE: the legacy-runner guard against a SELECTOR-provided remote executor lives in
    // AgentResumeIntegrationTest, next to its plain-scope sibling.
}
