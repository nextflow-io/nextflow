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

package nextflow.processor

import com.google.common.hash.HashCode
import nextflow.Session
import org.pf4j.ExtensionPoint

/**
 * Plugin extension point deciding how a task is resolved against the cache: resumed from a previous
 * execution, or given a work directory and launched.
 *
 * <p>{@link TaskProcessor} resolves the registered strategies once, in {@link nextflow.plugin.Priority}
 * order, and uses the first one that {@link #isEnabled applies} to the session; with none registered,
 * or none applying, the {@link DefaultTaskCacheStrategy} is used, so a run without such a plugin
 * resolves its tasks exactly as before.
 *
 * <p>A strategy is expected to abstain (return {@code false} from {@link #isEnabled}) whenever its cache
 * is not the one in use for the session: the plugin registry is process-wide, while the choice of
 * strategy belongs to the run.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
interface TaskCacheStrategy extends ExtensionPoint {

    /**
     * Whether this strategy applies to the run; the highest-priority applicable one wins, else the
     * default.
     *
     * @param session The session of the run.
     * @return {@code true} to take over the resolution of every task of the run.
     */
    boolean isEnabled(Session session)

    /**
     * Resolve {@code task}, whose hash is {@code hash}: end by calling exactly one of
     * {@link TaskResolver#resume resolver.resume} (returning {@code true}) or
     * {@link TaskResolver#launch resolver.launch}.
     *
     * @param task The task to resolve.
     * @param hash The task hash, as computed by its {@link TaskHasher}.
     * @param tryCache Whether a cached execution may be resumed; {@code false} on a retry, and when
     *      the run is not resuming.
     * @param resolver The processor's primitives to resume or launch the task with.
     */
    void resolve(TaskRun task, HashCode hash, boolean tryCache, TaskResolver resolver)

}
