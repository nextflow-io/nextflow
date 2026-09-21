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
import groovy.transform.Immutable

/**
 * Resolved identity of an agent task, attached to the synthetic process config by
 * {@code AgentDef.buildAgentTask} so an observer can tell an agent task from an ordinary
 * process task and record what the agent actually was.
 *
 * <p>The carrier is a plain immutable object, NOT a Map or Closure, on purpose:
 * {@code TaskConfig} extends {@code LazyMap}, whose {@code get()} deep-copies Map values
 * and <em>invokes</em> Closure values on every access (and a Closure value would flip
 * {@code LazyMap.dynamic}, short-circuiting {@code TaskRun.hasCacheableValues()}). An
 * immutable POJO is returned by identity and is inert.
 *
 * <p>This carries only the values resolved at task-build time. The two runtime-only facts
 * live elsewhere: the concrete model reported by the provider travels through
 * {@link AgentCallInfo} into the task context key {@code $agentResolvedModel}, and the
 * rendered prompt is a body-closure local that is deliberately never persisted.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Immutable
@CompileStatic
class AgentTaskInfo {

    /**
     * Process-config key under which the info object is attached. Read back with an
     * {@code instanceof} guard, never a truthiness check, so that a user-supplied
     * {@code agentInfo} directive in a {@code withName:} selector is inert rather than a
     * way to forge an agent lineage record.
     */
    public static final String CONFIG_KEY = 'agentInfo'

    /** Name of the {@code AgentRunner} implementation selected for the run */
    String runner
    /** Effective model id, after the `agent.defaultModel` config fallback */
    String model
    /** Resolved `instruction:` directive */
    String instruction
    /** Resolved `goal:` directive */
    String goal
    /** Source text of the `prompt:` template (NOT the per-task rendered prompt) */
    String promptTemplate
    /** Effective iteration ceiling, after the `agent.maxIterationsDefault` fallback */
    int maxIterations
    /** Key-sorted JSON of the synthesized output schema; null for a free-text agent */
    String outputSchema
    /** Names of the resolved tools the agent may call; null when the agent declares none */
    List<String> tools
    /** Names of the resolved skills available to the agent; null when the agent declares none */
    List<String> skills
}
