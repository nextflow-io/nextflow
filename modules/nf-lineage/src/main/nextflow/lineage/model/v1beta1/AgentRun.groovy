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

package nextflow.lineage.model.v1beta1

import groovy.transform.Canonical
import groovy.transform.CompileStatic
import nextflow.lineage.serde.LinSerializable

/**
 * Models an agent execution: the sibling of {@link TaskRun} for the {@code agent} primitive.
 *
 * <p>An agent lowers to an ordinary task, so this record is written in place of a
 * {@link TaskRun} under the same task hash, and its results are recorded as a
 * {@link TaskOutput} exactly like a process. What differs is the identity being captured:
 * a process is identified by its script, an agent by the model it called and the tools,
 * skills and prompt template it called it with.
 *
 * <p>Deliberately NOT recorded, and why:
 * <ul>
 *   <li><b>the rendered prompt</b> — a body-closure local that is never persisted; it is
 *     derivable from {@code promptTemplate} plus the recorded {@code input}, and recording it
 *     would put interpolated user data (potentially secrets) on disk uncapped.</li>
 *   <li><b>the resolved script</b> — an agent has none worth keeping, and on the RPC path the
 *     rendered command embeds a per-invocation capability token that must not be persisted.</li>
 *   <li><b>token usage, turn count, tool-call count</b> — not instrumented anywhere in the
 *     runtime today; there is nothing to read at task-complete time.</li>
 * </ul>
 * Add fields at the END: {@code @Canonical} makes the declaration order the positional
 * constructor signature.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Canonical
@CompileStatic
class AgentRun implements LinSerializable {
    /**
     * Execution session identifier
     */
    String sessionId
    /**
     * Agent task name
     */
    String name
    /**
     * Checksum of the canonical agent identity source, i.e. the same text that feeds the
     * resume cache key: runner, model, instruction, goal, iteration ceiling, prompt template,
     * output schema and skills fingerprint
     */
    Checksum codeChecksum
    /**
     * Name of the agent runner implementation that executed the run
     */
    String runner
    /**
     * Model requested for the run, after the `agent.defaultModel` config fallback
     */
    String model
    /**
     * Concrete model reported by the provider, when available. May differ from {@link #model}
     * when a floating alias resolves to a dated snapshot. Null for a tool agent (which is
     * non-cacheable and so stores no context) and on the RPC runner path.
     */
    String resolvedModel
    /**
     * Resolved `instruction:` directive
     */
    String instruction
    /**
     * Resolved `goal:` directive
     */
    String goal
    /**
     * Source text of the `prompt:` template. Combined with {@link #input} this reconstructs
     * what the model was asked.
     */
    String promptTemplate
    /**
     * Effective iteration ceiling for the agent loop
     */
    int maxIterations
    /**
     * Key-sorted JSON of the synthesized output schema. Stored as text rather than a Map so a
     * round-trip cannot coerce schema integers to floating point.
     */
    String outputSchema
    /**
     * Names of the tools the agent was allowed to call
     */
    List<String> tools
    /**
     * Names of the skills made available to the agent
     */
    List<String> skills
    /**
     * Agent run input
     */
    List<Parameter> input
    /**
     * Container used for the agent run, when the runner executes out of process
     */
    String container
    /**
     * Workflow run associated to the agent run
     */
    String workflowRun
    /**
     * Remote Nextflow module that defines the agent executed by this run, encoded as
     * {@code name@version}. Null when the agent is not defined in a remote module.
     */
    String moduleId
}
