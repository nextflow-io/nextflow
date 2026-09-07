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
 * Core-owned seam for carrying the concrete model snapshot resolved by an
 * {@link AgentRunner} back to the core task body, without violating the
 * core↔plugin boundary (the plugin still returns only a {@code String} from
 * {@code run}). Mirrors the {@code ModuleToolBridge.CONTEXT} ThreadLocal
 * pattern: the plugin writes via {@link #setResolvedModel} on the same
 * exec-service body thread that then reads via {@link #consumeResolvedModel}.
 *
 * <p>Used for resume drift observability (design §9.5/D6): after
 * {@code model.chat(...)} the plugin stashes {@code response.metadata().modelName()};
 * the core body reads it once and stores it in the task context under
 * {@code $agentResolvedModel} so it is persisted/replayed by the cache.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class AgentCallInfo {

    private static final ThreadLocal<String> RESOLVED_MODEL = new ThreadLocal<String>()

    /** Set the concrete model snapshot resolved on the current thread. Called by the plugin after {@code model.chat}. */
    static void setResolvedModel(String modelName) {
        RESOLVED_MODEL.set(modelName)
    }

    /**
     * Read and clear the resolved model snapshot for the current thread.
     * Returns {@code null} when none was set (e.g. no fresh call happened).
     */
    static String consumeResolvedModel() {
        final value = RESOLVED_MODEL.get()
        RESOLVED_MODEL.remove()
        return value
    }

    /** Clear any resolved model snapshot for the current thread. */
    static void clear() {
        RESOLVED_MODEL.remove()
    }
}
