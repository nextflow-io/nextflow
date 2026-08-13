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

package nextflow.agent.rpc

import groovy.transform.CompileStatic
import nextflow.Session
import nextflow.SysEnv
import nextflow.config.spec.ConfigOption
import nextflow.config.spec.ConfigScope
import nextflow.script.dsl.Description
import nextflow.util.Duration
import nextflow.agent.AgentConfig
import nextflow.agent.AgentRunnerRequest

/**
 * Model the nested `agent.rpc` scope: the driver-side broker a canonical agent task dials back to.
 *
 * <p>Declared as a nested {@link ConfigScope} rather than a plain map option so each key is
 * individually validated -- see {@link AgentConfig#AGENT_ONLY_OPTIONS}. Read once per session by
 * the runner plugin's {@code AgentRpcBroker} (no longer a core class), so a value written inside a
 * selector block resolves but has no effect.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class AgentRpcConfig implements ConfigScope {

    static final int DEFAULT_PORT = 0

    /**
     * The name a container gets, per container ENGINE, for the host its driver runs on. Only these
     * two engines synthesize such a name, which is why the alias is one rung of the ladder and not
     * the whole of it: every other engine, and every executor that launches the task somewhere other
     * than the driver host, is answered by inference instead -- see
     * {@link AgentRpcHostResolver}, which owns the ladder, and which rejects the cases that cannot
     * work before the run starts rather than advertising an address the task cannot reach.
     *
     * <p>The alias is also WRONG, not merely absent, when the driver is ITSELF containerized: it
     * names the host, while the broker listens inside the driver container. That case is row R3 of
     * the ladder and is why the engine name alone is not enough to answer this question.
     *
     * <p>NOTE for {@code docker}: the name is built in on Docker Desktop, while Linux Docker
     * (>= 20.10) resolves it only when the container is run with
     * {@code --add-host=host.docker.internal:host-gateway};
     * {@link nextflow.agent.AgentLaunchConditions#withDockerHostGateway} adds that run option
     * to the agent task so both behave the same.
     */
    private static final Map<String,String> ENGINE_HOST_ALIASES = Collections.unmodifiableMap(
        [ docker: 'host.docker.internal', podman: 'host.containers.internal' ] as Map<String,String> )

    /** The built-in host alias for a container engine, or {@code null} when the engine has none. */
    static String hostAliasFor(String engine) {
        return engine != null ? ENGINE_HOST_ALIASES.get(engine) : null
    }

    @ConfigOption
    @Description("""
        The port the driver's agent RPC broker binds to. `0` (the default) picks an ephemeral port; pin it when the driver sits behind a firewall that must be opened in advance.
    """)
    final Integer port

    /**
     * Environment fallback for {@link #remoteHost}. The driver's address is a property of WHERE the
     * pipeline runs, not of the pipeline, so it belongs in the deployment's environment as readily
     * as in its config: the same script can then run on a laptop and in a cluster without editing
     * a config file that is otherwise portable.
     */
    static final String REMOTE_HOST_ENV = 'NXF_AGENT_RPC_REMOTE_HOST'

    @ConfigOption
    @Description("""
        The host name a containerized agent task uses to reach the driver's RPC broker. Falls back to the `NXF_AGENT_RPC_REMOTE_HOST` environment variable, then to inference from the container engine and the executor: the engine's built-in host alias (`host.docker.internal` for Docker, `host.containers.internal` for Podman) when that engine runs on the driver host, `127.0.0.1` when the container shares the driver's network namespace, and otherwise the driver's own address on its default route. Set it explicitly where the driver does not share a network with its tasks: a Kubernetes driver outside the cluster the pods run in, a cloud batch driver on an instance in a different VPC from the compute environment, or a multi-homed submit node whose default route is not the compute fabric. Give a routable driver address or an in-cluster service name.
    """)
    final String remoteHost

    /**
     * Which of the two explicit rungs answered {@link #remoteHost} -- {@code agent.rpc.remoteHost}
     * or {@code NXF_AGENT_RPC_REMOTE_HOST} -- or {@code null} when neither did and the address is
     * therefore inferred. Kept because the source is what the registration line reports (see
     * {@link AgentRpcHost#describe}): the expensive failure of an inferred address is that it looks
     * exactly like a configured one, so an operator must be able to see which they got.
     */
    final String remoteHostSource

    /**
     * How long an unconsumed agent RPC capability stays valid. This clock starts when the task
     * SCRIPT is generated, not when the job runs, so it has to absorb the executor's queueing
     * latency: one hour covers essentially any real queue wait, while a capability remains
     * single-use and unguessable for that window.
     */
    static final Duration DEFAULT_CAPABILITY_TIMEOUT = Duration.of('1h')

    @ConfigOption
    @Description("""
        How long an agent RPC capability remains valid while its task waits to start, i.e. the maximum queueing delay tolerated between the task script being generated and the agent connecting back to the driver (default: `1 h`).
    """)
    final Duration capabilityTimeout

    @ConfigOption
    @Description("""
        When `true` (the default), the driver's agent RPC broker serves TLS with a certificate generated for the run, whose fingerprint the agent task pins. Set it to `false` only to inspect the stream while debugging: the prompt, inputs, tool arguments and results then cross the network in cleartext, and the task cannot tell the real driver from an impostor.
    """)
    final Boolean tls

    /* required by the spec reflection -- do not remove */
    AgentRpcConfig() {}

    AgentRpcConfig(Map opts) {
        port = opts.port != null ? opts.port as Integer : DEFAULT_PORT
        // NOT defaulted here: the remaining rungs need the container engine, the executor and facts
        // about the driver host that this object cannot see -- see AgentRpcHostResolver
        remoteHost = resolveConfiguredHost(opts, SysEnv.get())
        remoteHostSource = remoteHost == null ? null
                : (opts?.remoteHost?.toString() ? 'agent.rpc.remoteHost' : REMOTE_HOST_ENV)
        capabilityTimeout = opts.capabilityTimeout != null ? opts.capabilityTimeout as Duration : DEFAULT_CAPABILITY_TIMEOUT
        tls = opts.tls != null ? opts.tls as Boolean : Boolean.TRUE
    }

    /**
     * The first two rungs of the ladder -- {@code agent.rpc.remoteHost}, then
     * {@code NXF_AGENT_RPC_REMOTE_HOST} -- resolved against an explicit environment so tests can
     * swap it instead of mutating the process. Config wins: a value written for this pipeline is
     * more specific than one exported for whatever else shares the shell.
     *
     * <p>Each rung is tested for TRUTHINESS, so an empty value (an exported-but-unset variable, or
     * {@code remoteHost = params.host} with {@code params.host} absent) falls through instead of
     * shadowing the rungs below it and advertising an empty host.
     */
    protected static String resolveConfiguredHost(Map opts, Map env) {
        final configured = opts?.remoteHost?.toString()
        if( configured )
            return configured
        final fromEnv = env?.get(REMOTE_HOST_ENV)?.toString()
        return fromEnv ?: null
    }

    Integer getPort() { port }

    String getRemoteHost() { remoteHost }

    String getRemoteHostSource() { remoteHostSource }

    Duration getCapabilityTimeout() { capabilityTimeout }

    Boolean getTls() { tls }

    /**
     * Whether the broker serves TLS; on unless explicitly disabled, so the no-arg constructor used
     * by the spec reflection also reports the secure default.
     */
    boolean tlsEnabled() { tls == null || tls }

    /**
     * The host name to advertise to a containerized agent task launched by the given engine ON THE
     * DRIVER HOST -- a THIN CALLER of {@link AgentRpcHostResolver}, which owns the ladder, so this
     * and {@link nextflow.agent.AgentLaunchConditions#requireBrokerHost} cannot drift apart.
     * {@code null} when no row applies, which the guard has already rejected before the run
     * started.
     *
     * <p>Only the LOCAL rows can be answered from an engine name alone, which is exactly what this
     * signature carries; the broker resolves through {@link #resolveBrokerHost} instead, because a
     * session-level engine name says nothing about an executor that runs the task elsewhere.
     */
    String resolveRemoteHost(String engine, Session session = null) {
        final result = AgentRpcHostResolver.of(session)
                .resolve(null, AgentConfig.DEFAULT_EXECUTOR, engine, session?.getContainerConfig(engine), null, this)
        return result.resolved ? result.host : null
    }

    /**
     * The address the broker falls back to when the registration carries none of its own.
     *
     * <p>The address that is actually advertised rides on {@link AgentRunnerRequest#brokerHost}: the
     * pre-ignition guard resolves it PER AGENT DEFINITION, with the full context (the executor
     * instance, the engine config, the task's container options), and a run may legitimately hold
     * several -- one agent on the local docker engine and another on {@code k8s} resolve different
     * addresses, and each task must be told its own. None of that context is recoverable here.
     *
     * <p>This is therefore only for a runner that registers WITHOUT the guard, or a spec: it answers
     * from what the session alone knows, which is exactly what shipped behaviour did -- the local
     * rows, keyed off the enabled engine. It can be an error row, and {@code register()} rejects with
     * that row's message rather than advertising {@code null:<port>}.
     */
    AgentRpcHost resolveBrokerHost(Session session) {
        final engine = session?.getContainerConfig()
        return AgentRpcHostResolver.of(session)
                .resolve(null, AgentConfig.DEFAULT_EXECUTOR, engine?.getEngine(), engine, null, this)
    }
}
