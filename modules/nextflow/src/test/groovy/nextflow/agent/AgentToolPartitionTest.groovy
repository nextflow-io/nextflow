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
import nextflow.processor.TaskConfig
import nextflow.script.AgentBuilder.AgentInput
import nextflow.script.AgentBuilder.AgentOutput
import nextflow.script.AgentDef
import nextflow.script.BaseScript
import nextflow.script.PromptDef
import nextflow.script.ScriptBinding
import spock.lang.Timeout
import test.Dsl2Spec
import test.MockSession
import nextflow.agent.rpc.AgentRpcRegistration

/**
 * The §5 runner split, end to end through the real agent-build path.
 *
 * <p>Every assertion here is about the PARTITION rather than about any single tool: the resolved
 * selection leaves {@code AgentDef} in two disjoint halves — the brokered descriptors on
 * {@code toolSpecs}, the runner-native wire names on {@code nativeToolNames} — and each half must
 * arrive intact at the runner that serves it. The two failure modes this pins down are both
 * SILENT, which is why they need their own tests: a native name that never reaches the request
 * leaves the model with no tool and no error, and a native name that reaches {@code toolSpecs}
 * enters the broker's allowlist and gets executed in the driver JVM instead of the container.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Timeout(30)
class AgentToolPartitionTest extends Dsl2Spec {

    def cleanup() {
        AgentRunnerProvider.testRunner = null
    }

    // -----------------------------------------------------------------------
    // §5 invariant — a containerized runner
    // -----------------------------------------------------------------------

    def 'the runner-native names reach a canonical runner on their own field, never as descriptors'() {
        given: 'a launch-spec runner, i.e. one that owns a container -- what `shell:` requires'
        final List<AgentRunnerRequest> registered = []
        AgentRunnerProvider.testRunner = canonicalRunner(registered)
        newSession(containerized())

        when: 'the agent selects a partial fs: set plus the shell family'
        runTaskBody(newAgent([model: 'openai/gpt-4o', tools: ['fs:read', 'fs:grep', 'shell:bash']]),
                new TaskConfig([container: 'agent-image:test']))

        then: 'the native half travels beside toolSpecs, in inventory order'
        registered.size() == 1
        registered[0].nativeToolNames == ['read', 'grep', 'bash']

        and: 'and NEVER inside it: no descriptor is minted for a tool the runner serves itself'
        !registered[0].toolSpecs
        registered[0].brokeredToolNames().isEmpty()

        and: 'the partition is enforced, not merely observed'
        registered[0].checkToolPartition()
    }

    def 'a canonical agent selecting only fs: gets no sandbox context in the driver'() {
        given: 'the container roots the fs: builtins at its own cwd -- the driver has no business there'
        final List<AgentRunnerRequest> registered = []
        AgentRunnerProvider.testRunner = canonicalRunner(registered)
        newSession(containerized())

        when:
        runTaskBody(newAgent([model: 'openai/gpt-4o', tools: ['fs:*']]),
                new TaskConfig([container: 'agent-image:test']))

        then: 'nothing is dispatchable back into the driver'
        registered[0].nativeToolNames == ['read', 'write', 'edit', 'ls', 'grep', 'find']
        !registered[0].toolSpecs
        and: 'a driver-side fs: call is impossible, so a remote work dir cannot break one'
        registered[0].dispatch == null || registered[0].dispatch.call('read', '{"path":"x"}').contains('Unknown agent tool')
    }

    // -----------------------------------------------------------------------
    // the in-JVM runner serves the SELECTED leaves, and only those
    // -----------------------------------------------------------------------

    def 'an in-JVM agent declaring one fs: leaf is served that leaf and no other'() {
        given:
        AgentRunnerRequest captured = null
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req -> captured = req; 'ok' } as AgentRunner
        newSession()

        when: 'read-only access is declared'
        runTaskBody(newAgent([model: 'openai/gpt-4o', tools: ['fs:read']]), new TaskConfig([:]))

        then: 'exactly one native name travels -- a partial selection is never widened to fs:*'
        captured.nativeToolNames == ['read']

        and: 'the declared leaf IS served in the driver JVM (it fails on the missing work dir, not on the name)'
        captured.dispatch.call('read', '{"path":"x"}').contains('`read` tool unavailable')

        and: 'while the leaves that were NOT declared are not tools at all -- no silent write access'
        captured.dispatch.call('write', '{"path":"x","content":"y"}').contains('Unknown agent tool `write`')
        captured.dispatch.call('edit', '{"path":"x"}').contains('Unknown agent tool `edit`')
    }

    // -----------------------------------------------------------------------
    // `shell:` is pi-only (§5), and the refusal must name the runner that has it
    // -----------------------------------------------------------------------

    def 'shell:bash is refused at agent-build time on an in-JVM runner, naming pi'() {
        given:
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req -> 'never reached' } as AgentRunner
        newSession()

        when:
        newAgent([model: 'openai/gpt-4o', tools: ['shell:bash']]).buildAgentTask(['hello'])

        then: 'the error names the ref, the reason and the runner that can serve it'
        def e = thrown(Exception)
        e.message.contains('`shell:bash`')
        e.message.contains('pi')
    }

    // -----------------------------------------------------------------------
    // helpers (mirrors AgentAsTaskIntegrationTest)
    // -----------------------------------------------------------------------

    /** The minimal config that containerizes a canonical agent task: an engine plus an image. */
    private static Map containerized(Map config = [:]) {
        final result = new LinkedHashMap(config)
        result.docker = [enabled: true]
        final agentScope = new LinkedHashMap((Map) (config.agent ?: [:]))
        agentScope.container = 'agent-image:test'
        result.agent = agentScope
        return result
    }

    /** A launch-spec runner recording every registered request, so the wire contract is assertable. */
    private static AgentRunner canonicalRunner(List<AgentRunnerRequest> registered) {
        return new AgentRunner() {
            @Override
            String getName() { 'external' }

            @Override
            AgentLaunchSpec getLaunchSpec() {
                new AgentLaunchSpec(
                    containerProxyCommand: ['/opt/agent-rpc'],
                    containerHarnessCommand: ['node', '/opt/runner.mjs'])
            }

            @Override
            AgentRpcRegistration register(AgentRunnerRequest request, boolean remote) {
                registered << request
                return new AgentRpcRegistration('inv-1', 'tok-1', 'host.docker.internal:9999', 'abc123')
            }

            @Override
            String run(AgentRunnerRequest request) { throw new UnsupportedOperationException('canonical task path') }
        }
    }

    private Session newSession(Map config = null) {
        def session = config ? new MockSession(config) : new MockSession()
        session.setBinding(new ScriptBinding())
        session.init(null, null, null, null)
        session.start()
        Global.session = session
        return session
    }

    /** Lower the agent and invoke the synthesized body against a stand-in task context. */
    private static String runTaskBody(AgentDef agent, TaskConfig taskConfig) {
        final body = (Closure) agent.buildAgentTask(['hello']).getTaskBody().closure.clone()
        body.setDelegate([q: 'hello', task: taskConfig])
        body.setResolveStrategy(Closure.DELEGATE_ONLY)
        return body.call()
    }

    private AgentDef newAgent(Map directives = [model: 'openai/gpt-4o']) {
        final owner = Mock(BaseScript) { getBinding() >> new ScriptBinding() }
        return new AgentDef(owner, 'qa', directives as Map<String,Object>,
            [new AgentInput('q', String)], [new AgentOutput('answer', String)],
            new PromptDef({ -> 'Q' }, 'Q'))
    }
}
