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

import java.nio.file.Path

import com.google.common.hash.HashCode

/**
 * The primitives a {@link TaskCacheStrategy} needs to resolve a task against the cache -- and
 * nothing else of the {@link TaskProcessor} that implements it.
 *
 * <p>A strategy decides <i>which</i> cache entry to look at and <i>which</i> work directory a task
 * executes in; the resolver carries out the two possible outcomes -- resume the task from an entry,
 * or launch it in a work directory -- with the processor's own machinery, so a strategy never
 * touches the processor, the session or the cache store directly.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
interface TaskResolver {

    /**
     * Look up the cache entry stored under the given key.
     *
     * @param key The task hash the entry is keyed by.
     * @return The entry, or {@code null} when the cache has none for {@code key}.
     */
    TaskEntry entry(HashCode key)

    /**
     * Resume a task from a cache entry: verify the outputs in {@code workDir} against the entry and,
     * when they are all there, bind them and announce the cache hit.
     *
     * <p>The resume works on a <b>copy</b> of {@code task}: the copy becomes the cached task (its work
     * dir, hash, context and outputs are set), while {@code task} itself is left untouched and stays
     * launchable, which every strategy relies on when a failed resume falls through to a launch.
     *
     * @param task The task instance to resume into.
     * @param key The hash the task is resumed under, i.e. the key its entry was found by.
     * @param workDir The work directory holding the cached outputs.
     * @param entry The cache entry to resume from.
     * @return {@code true} when the task was resumed, {@code false} when the entry could not be used
     *      (missing outputs, non-zero exit) and the task still has to be launched.
     */
    boolean resume(TaskRun task, HashCode key, Path workDir, TaskEntry entry)

    /**
     * Launch a task for execution under the given hash in the given work directory.
     *
     * @param task The task to execute.
     * @param key The hash the task executes under; it names the work directory and keys the entry
     *      written on completion.
     * @param workDir The work directory the task executes in.
     */
    void launch(TaskRun task, HashCode key, Path workDir)

    /**
     * The work directory a hash maps to under the executor's work root -- where a task executing under
     * {@code key} would run.
     *
     * @param key A task hash.
     * @return The work directory of {@code key}; nothing is created.
     */
    Path workDirFor(HashCode key)

}
