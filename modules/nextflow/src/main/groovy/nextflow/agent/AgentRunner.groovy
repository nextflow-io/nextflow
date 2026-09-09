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
import nextflow.agent.rpc.AgentRpcHostResolver
import nextflow.agent.rpc.AgentRpcRegistration
import nextflow.agent.rpc.AgentRpcConfig

/**
 * SPI implemented by an agent runner plugin (e.g. nf-agent). Given a resolved
 * {@link AgentRunnerRequest}, drive the LLM and return the final assistant text.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
interface AgentRunner {

    /**
     * Stable user-facing runner identifier. Implementations supplied by plugins
     * should override this value (for example {@code pi} or {@code langchain4j}).
     * The default keeps this interface single-abstract-method compatible with
     * closure-coerced runners used by tests and embedding applications.
     */
    default String getName() { getClass().getSimpleName() }

    /**
     * Return a canonical task launch description, or {@code null} for a legacy
     * in-JVM runner. The default preserves closure-coerced test runners.
     */
    default AgentLaunchSpec getLaunchSpec() { null }

    /**
     * The container image a canonical task of this runner runs in when the user declares no
     * {@code agent.container}, or {@code null} when the runner has no image of its own.
     *
     * A runner that ships its runtime IN an image (see nf-agent-pi) generates this from its own
     * VERSION at build time, so the jar and the image it needs cannot drift: the coordinate the
     * jar asks for is by construction the tag the release publishes.
     *
     * The default keeps this interface single-abstract-method compatible with closure-coerced
     * runners used by tests and embedding applications, and leaves a runner that has no image
     * of its own -- an in-JVM one, for instance -- requiring an explicit {@code agent.container}.
     */
    default String getDefaultContainer() { null }

    /**
     * Issue the connection material a canonical agent task needs to call back into
     * the driver, i.e. the broker endpoint plus a single-use capability token.
     *
     * A runner that returns a {@link #getLaunchSpec()} MUST implement this, because
     * the launch command is built from the returned registration. The broker
     * implementation and its transport dependencies live with the runner plugin so
     * that a Nextflow distribution carrying no agent runner carries no RPC stack.
     *
     * <p>CONTRACT for {@code remote=true}, which is how a canonical (launch-spec) task is
     * always registered: the HOST part of the returned endpoint MUST be
     * {@link AgentRunnerRequest#brokerHost} -- the address {@link AgentRpcHostResolver}
     * resolved for THIS agent definition -- falling back to
     * {@link AgentRpcConfig#resolveBrokerHost} only when the request carries none. It is
     * per definition and not per run because a script may declare one agent on the local
     * docker engine and another on {@code k8s}, which resolve different addresses. Core
     * resolves and validates it on the runner's behalf before the run starts -- see
     * {@link nextflow.agent.AgentLaunchConditions#requireBrokerHost} -- so a runner that
     * derives the driver address by some other means would advertise one address while the
     * run was admitted on another, and the cases core rejects pre-ignition would fail an
     * hour later instead, when the capability finally expires.
     *
     * <p>CONTRACT on transport security: the returned registration MUST carry either the SHA-256
     * fingerprint of the broker's certificate, or {@code insecure = true} to declare that the broker
     * deliberately serves cleartext. There is no third state. A registration that carries neither is
     * rejected by {@link AgentRpcRegistration#transportArgs()} at SCRIPT-GENERATION time, before any
     * job is submitted -- deliberately, because the alternative is a proxy told to dial cleartext at
     * a TLS listener, which surfaces inside the task as an unrelated connection failure. A runner
     * that predates transport security, and returns only invocation id, token and endpoint, fails
     * here and must be updated rather than configured around.
     *
     * The default keeps this interface single-abstract-method compatible with
     * closure-coerced runners used by tests and embedding applications.
     */
    default AgentRpcRegistration register(AgentRunnerRequest request, boolean remote) {
        throw new UnsupportedOperationException("Agent runner `${getName()}` provides a launch spec but no RPC broker")
    }

    String run(AgentRunnerRequest request)
}
