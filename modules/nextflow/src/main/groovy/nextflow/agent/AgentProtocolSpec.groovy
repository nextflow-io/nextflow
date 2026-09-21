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

/**
 * Build the portable request payload shared by agent protocol transports.
 *
 * <p>INVARIANT: this payload carries NO credential. It is the PORTABLE half of a start frame --
 * relayed verbatim by the agent proxy, and the shape a transport is free to log or persist -- so
 * the credential resolved for the request travels BESIDE it, as a top-level frame field the RPC
 * broker adds only when the link is TLS-protected (see {@code AgentRpcBroker#startFrame}), and
 * {@link AgentRunnerRequest#apiKey} stays deliberately omitted here. {@code baseUrl} is not a
 * secret and does travel, because the remote runner must target the same endpoint the driver
 * resolved.
 *
 * <p>The agent's tools cross in TWO fields, and the split is load-bearing: {@code toolSpecs}
 * carries the brokered descriptors the runner calls back into the driver for, while
 * {@code nativeToolNames} carries bare names the runner is told to enable from its OWN tool set
 * (on pi, the SDK builtins added to the session allowlist). This is the single point at which both
 * halves are handed to a runner, so it is where the partition between them is checked.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
final class AgentProtocolSpec {

    private AgentProtocolSpec() {}

    static Map<String,Object> fromRequest(AgentRunnerRequest request) {
        // A native name that also appears among the brokered descriptors would reach the runner as
        // a tool it must call the driver for -- the exact confusion the split exists to prevent.
        request.checkToolPartition()
        return [
            model: request.model,
            instruction: request.instruction,
            goal: request.goal,
            prompt: request.prompt,
            inputJson: request.inputJson,
            outputSchema: request.outputSchema,
            toolSpecs: request.toolSpecs,
            nativeToolNames: request.nativeToolNames,
            skills: request.skills,
            maxIterations: request.maxIterations > 0 ? request.maxIterations : 20,
            trace: request.trace,
            temperature: request.temperature,
            workDir: request.workDir,
            // NOTE: `apiKey` is intentionally NOT part of this payload -- see the class javadoc
            baseUrl: request.baseUrl ]
    }
}
