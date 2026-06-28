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
 * Being {@code @Canonical}, the positional constructor order follows the field
 * declaration order below:
 * {@code (model, instruction, prompt, maxIterations, tools, outputSchema, inputJson, toolSpecs, dispatch, requestTimeoutSeconds, goal, agentName, trace, skills)}.
 * In practice the request is built with named arguments, so the order is not
 * relied upon at the call site.
 *
 * The {@code requestTimeoutSeconds} carries the configured per-request LLM chat
 * timeout (from the {@code agent.requestTimeout} config option); when {@code 0}
 * the runner applies its own built-in default.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Canonical
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
    /** Portable descriptors of the agent's declared skills (langchain4j-free); the
     * runner maps them onto langchain4j {@code Skill}s and exposes them via the
     * skills tool provider. Empty/null when no {@code skills} directive is declared. */
    List<SkillDescriptor> skills
}
